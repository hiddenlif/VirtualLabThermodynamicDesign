"""Whole-domain thermodynamic enclosure equations.

Dyadic trees cover the declared parameter box by construction. An exclusion
leaf names an observation whose interval cannot intersect the observation band.
A retained leaf includes an operating-domain tree proving the lower bound.
All splits are reconstructed by the checker. Failed proof construction is
INCONCLUSIVE, never scientific insufficient evidence or rejection.
"""
from __future__ import annotations
import math
import time
from copy import deepcopy

from mpmath import iv

iv.dps = 40
R = iv.mpf(8.31446261815324)
T0 = iv.mpf(298.15)
LIGAND = iv.mpf(2e-6)


def interval(lo, hi=None):
    return iv.mpf([lo, lo if hi is None else hi])


def endpoints(x):
    # Outward binary64 conversion of arbitrary precision interval endpoints.
    return (math.nextafter(float(x.a), -math.inf),
            math.nextafter(float(x.b), math.inf))


def point_logk(log0, h0, cp, t):
    return log0 + h0/R*(1/T0-1/t) + cp/R*(iv.ln(t/T0)+T0/t-1)


def logk(log0, h0, cp, temp):
    # Tight enclosure of the scalar temperature basis; B has minimum at T0.
    tl, th = endpoints(temp)
    a = 1/T0 - 1/temp
    def b(t):
        return iv.ln(t/T0) + T0/t - 1
    b1, b2 = b(interval(tl)), b(interval(th))
    bl = 0.0 if tl <= 298.15 <= th else min(endpoints(b1)[0], endpoints(b2)[0])
    bh = max(endpoints(b1)[1], endpoints(b2)[1])
    return log0 + h0/R*a + cp/R*interval(bl, bh)


def thermo(box, ph_range, temp_range, buffer_h):
    p = [interval(*b) for b in box]
    ph, t = interval(*ph_range), interval(*temp_range)
    ln10 = iv.ln(interval(10))
    lf = logk(ln10*p[3], interval(-25000), interval(-100), t) - ln10*ph
    lc = logk(ln10*p[4], interval(-32000), interval(-140), t) - ln10*ph
    uf, uc = iv.exp(lf), iv.exp(lc)
    def fraction(u):
        a,b=endpoints(u);low=interval(a)/(1+interval(a));high=interval(b)/(1+interval(b))
        return interval(endpoints(low)[0],endpoints(high)[1])
    ff, fc = fraction(uf), fraction(uc)
    k = iv.exp(logk(p[0], p[1], p[2], t)) * (1+uc)/(1+uf)
    h = p[1]+p[2]*(t-T0)+fc*(interval(buffer_h)-32000-140*(t-T0))-ff*(interval(buffer_h)-25000-100*(t-T0))
    return k, h


def root_point(a, b, k):
    if b == 0 or a == 0:
        return interval(0)
    a, b, k = interval(a), interval(b), interval(k)
    s = a+b+1/k
    # Nonnegative expression avoids cancellation in the discriminant.
    discriminant = (a-b)**2 + 2*(a+b)/k + 1/k**2
    return 2*a*b/(s+iv.sqrt(discriminant))


def complex_box(a, b, k):
    # The physical root is monotone in both positive totals and association K.
    al, ah = endpoints(a); bl, bh = endpoints(b); kl, kh = endpoints(k)
    low = root_point(max(0, al), max(0, bl), kl)
    high = root_point(ah, bh, kh)
    return interval(max(0, endpoints(low)[0]), endpoints(high)[1])


def predict_box(box, action):
    """Displacement-corrected cumulative heats, including the shared offset."""
    if action['kind']=='derived':
        a=action['source'];k,h=thermo(box,[a['ph']]*2,[a['temperature']]*2,a['buffer_h']);values=[iv.ln(k),h];return [(v-action['centers'][i])/action['scales'][i] for i,v in enumerate(values)]
    if action['kind']=='blank':
        gain=interval(0);out=[]
        for inj in action['injections_l']:
            gain=gain*(1-interval(inj)/interval(action['volume_l']))+1
            out.append(interval(*box[7])*gain)
        return out
    if action['kind']=='standard':
        k=iv.exp(logk(iv.ln(interval(1e9)),interval(-35000),interval(0),interval(action['temperature'])))
        h=interval(-35000)
    else:k,h=thermo(box,[action['ph']]*2,[action['temperature']]*2,action['buffer_h'])
    v=interval(action['volume_l']);protein=interval(*box[5+action['lot']])*interval(action['protein_m']);ligand=interval(0);gain=interval(0);out=[]
    for inj in action['injections_l']:
        f=interval(inj)/v;protein*=1-f;ligand=ligand*(1-f)+f*interval(action['syringe_m']);gain=gain*(1-f)+1
        out.append(v*complex_box(protein,ligand,k)*h+interval(*box[7])*gain)
    return out


def response_box(box, operating, lot=0):
    k, _ = thermo(box, operating[0], operating[1], 0)
    kl = k*LIGAND
    # Monotone rational function, avoiding self-dependency expansion.
    low, high = endpoints(kl)
    lo = interval(low)/(1+interval(low))
    hi = interval(high)/(1+interval(high))
    return interval(*box[5+lot])*interval(endpoints(lo)[0], endpoints(hi)[1])


def predict_raw_box(box,action):
    if action['kind']=='derived':
        a=action['source'];k,h=thermo(box,[a['ph']]*2,[a['temperature']]*2,a['buffer_h']);values=[iv.ln(k),h];return [(v-action['centers'][i])/action['scales'][i] for i,v in enumerate(values)]
    if action['kind']=='blank':return [interval(*box[7]) for _ in action['injections_l']]
    if action['kind']=='standard':
        k=iv.exp(logk(iv.ln(interval(1e9)),interval(-35000),interval(0),interval(action['temperature'])));h=interval(-35000)
    else:k,h=thermo(box,[action['ph']]*2,[action['temperature']]*2,action['buffer_h'])
    v=interval(action['volume_l']);protein=interval(*box[5+action['lot']])*interval(action['protein_m']);ligand=old=interval(0);out=[]
    for inj in action['injections_l']:
        f=interval(inj)/v;protein*=1-f;ligand=ligand*(1-f)+f*interval(action['syringe_m']);new=complex_box(protein,ligand,k);out.append(v*(new-(1-f)*old)*h+interval(*box[7]));old=new
    return out

def decreasing_response(box,op):
    p=[interval(*b) for b in box];t=interval(*op[1]);ln10=iv.ln(interval(10))
    difference=logk(ln10*(p[4]-p[3]),interval(-7000),interval(-40),t)
    _,h=thermo(box,op[0],op[1],0)
    return endpoints(difference)[0]>=0 and endpoints(h)[1]<=0
