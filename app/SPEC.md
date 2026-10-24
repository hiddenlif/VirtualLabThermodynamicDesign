# Thermodynamic Calorimetry Experiment Design Specification

## 1. IPC Protocol & Execution Constraints
- The planner executes via `/app/run.sh` and communicates with the instrument over standard input and standard output using newline-delimited JSON.
- The instrument provides initial uncalibrated thermograms and campaign parameters.
- The planner may submit up to four sequential experimental batches (`{"op": "batch", "experiments": [...]}`).
- When planning is complete, the planner emits `{"op": "release"}` and writes the final certificate to `/app/release.json`.
- Per-campaign execution budget: 1800 seconds.

## 2. Parameter Coordinates
The physical model operates over eight continuous parameter coordinates:
1. `log_Kint`: Intrinsic association log-constant.
2. `H_int`: Intrinsic enthalpy of binding (J/mol).
3. `Cp_int`: Intrinsic heat capacity change (J/(mol*K)).
4. `pKa_free`: Free protein ionization pKa.
5. `pKa_bound`: Bound complex ionization pKa.
6. `active_0`: Active protein fraction for lot 0 (in [0, 1]).
7. `active_1`: Active protein fraction for lot 1 (in [0, 1]).
8. `heat_offset`: Instrument baseline heat offset.

## 3. Release Certificate Schema
The certificate written to `/app/release.json` must be a JSON object with:
- `witness`: List of 8 finite floats representing a parameter vector inside the prior domain that satisfies all acquired observation bands.
- `brackets`: A pair of 2-element intervals `[[lo_0, hi_0], [lo_1, hi_1]]` enclosing the worst-case operating response for each lot.
- `determinations`: A pair of string labels (`"CERTIFIED_COMPLIANT"` or `"CERTIFIED_VIOLATING"`) for each lot.
- `tree`: A verified domain-decomposition tree covering the entire eight-coordinate prior.

## 4. Certification & Tolerance Thresholds
- `claim_floor`: Target operating binding response threshold (0.52).
- `max_bracket_width`: Maximum allowed response bracket uncertainty (hi - lo <= 0.025).
- A lot is `"CERTIFIED_COMPLIANT"` only if `lo >= claim_floor`.
- A lot is `"CERTIFIED_VIOLATING"` only if `hi < claim_floor`.
- Tree complexity bounds: maximum node count and maximum subdivision depth specified in campaign configuration.

