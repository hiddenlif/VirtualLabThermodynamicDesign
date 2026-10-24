"""Analytic calibration and single-curve constraints, checked at endpoints."""

import copy
import math

import enclosures as E


def band(y, eps):
    return E.interval(y) + E.interval(-eps, eps)


def totals(action, active):
    p = E.interval(*active) * E.interval(action["protein_m"])
    l = E.interval(0)
    records = []
    for inj in action["injections_l"]:
        f = E.interval(inj) / E.interval(action["volume_l"])
        p *= 1 - f
        l = l * (1 - f) + f * E.interval(action["syringe_m"])
        records.append((p, l))
    return records


def corrected(history, box, j):
    h = history[j]
    s = E.interval(0)
    w = E.interval(0)
    out = []
    for y, inj in zip(h["y"], h["action"]["injections_l"]):
        f = E.interval(inj) / E.interval(h["action"]["volume_l"])
        s = s * (1 - f) + band(y, h["epsilon"])
        w = w * (1 - f) + 1
        out.append(s - E.interval(*box[7]) * w)
    return out


def ratio_at(box, action, logk):
    rs = totals(action, box[5 + action["lot"]])
    a, b = rs[0]
    c, d = rs[-1]
    k = E.iv.exp(E.interval(logk))
    # d log(C_first/C_last)/d log K = 1/S_first - 1/S_last.
    # S^2 = 1 + 2 K(P+L) + K^2(P-L)^2; these sufficient conditions
    # establish monotonic increase for every K>0 and every active fraction.
    if (
        E.endpoints(c + d - a - b)[0] <= 0
        or E.endpoints((c - d) ** 2 - (a - b) ** 2)[0] <= 0
    ):
        raise ValueError("Curve does not admit monotone ratio certificate")
    return E.complex_box(a, b, k) / E.complex_box(c, d, k)


def derive(box, history, j, log_range):
    h = history[j]
    a = h["action"]
    if a["kind"] != "binding":
        raise ValueError("Only binding curves have derived affinity constraints")
    q = corrected(history, box, j)
    if E.endpoints(q[0])[1] >= 0 or E.endpoints(q[-1])[1] >= 0:
        raise ValueError("Corrected heats must be strictly negative")
    ratio = q[0] / q[-1]
    low = ratio_at(box, a, log_range[0])
    high = ratio_at(box, a, log_range[1])
    if not low.b < ratio.a or not high.a > ratio.b:
        raise ValueError("Affinity endpoints do not exclude outer ranges")
    p, l = totals(a, box[5 + a["lot"]])[0]
    k = E.iv.exp(E.interval(*log_range))
    c = E.complex_box(p, l, k)
    H = q[0] / (E.interval(a["volume_l"]) * c)
    hr = E.endpoints(H)
    m = [sum(log_range) / 2, sum(hr) / 2]
    w = [max(m[0] - log_range[0], log_range[1] - m[0]), max(m[1] - hr[0], hr[1] - m[1])]
    w = [math.nextafter(math.nextafter(x, math.inf), math.inf) for x in w]
    # Virtual residuals are normalized independent bands; they do not reveal new data.
    return {
        "action": {
            "kind": "derived",
            "source": copy.deepcopy(a),
            "centers": m,
            "scales": w,
            "volume_l": 1.0,
            "injections_l": [1.0, 1.0],
        },
        "y": [0.0, 0.0],
        "epsilon": 1.0,
    }


def validate_calibration(old, history, new, blank_ids, standard_ids):
    if len(new) != 8 or any(new[i] != old[i] for i in range(5)):
        raise ValueError("Calibration only contracts active fractions and offset")
    for (a, b), (c, d) in zip(new, old):
        if not c <= a < b <= d:
            raise ValueError("Calibration is not an inner box")
    lo, hi = old[7]
    if not blank_ids:
        raise ValueError("Blank required")
    for j in blank_ids:
        h = history[j]
        if h["action"]["kind"] != "blank":
            raise ValueError("Calibration blank kind")
        for y in h["y"]:
            a, b = E.endpoints(band(y, h["epsilon"]))
            lo = max(lo, a)
            hi = min(hi, b)
    if not new[7][0] <= lo <= hi <= new[7][1]:
        raise ValueError("Offset interval omits possibilities")
    if len(standard_ids) != 2:
        raise ValueError("Two lot standards required")
    for lot, j in enumerate(standard_ids):
        h = history[j]
        if h["action"]["kind"] != "standard" or h["action"]["lot"] != lot:
            raise ValueError("Standard lot mismatch")
        target = band(h["y"][0], h["epsilon"])
        b = copy.deepcopy(new)
        b[5 + lot] = [new[5 + lot][0]] * 2
        lower = E.predict_raw_box(b, h["action"])[0]
        b[5 + lot] = [new[5 + lot][1]] * 2
        upper = E.predict_raw_box(b, h["action"])[0]
        # Known standard enthalpy is negative; cumulative heat decreases in r.
        if not lower.a > target.b or not upper.b < target.a:
            raise ValueError("Active-fraction endpoints fail monotone exclusion")
    return True
