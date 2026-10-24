Your task is to implement an autonomous calorimetry planner in `/app/run.sh`. The planner communicates with the instrument over standard input and standard output using newline-delimited JSON messages.

The instrument provides uncalibrated raw thermograms and campaign configuration parameters. The planner submits sequential titration batches, tracks irreversible protein and ligand stock depletion, and writes the final evaluation certificate to `/app/release.json`. Every lot response must be bounded by narrow bracket intervals, certified against the claim floor, and supported by a verified decomposition tree containing a compatible parameter witness.

Validate communication protocol handling and local execution by running `python3 /app/campaign.py` against `/app/practice_case.json`. Local runs test the protocol.

The forward linked-protonation model is implemented in `/app/physics.py`, and tree validation replay contracts reside in `/app/checker.py`, while parameter coordinate definitions, response bracket widths, and numerical verification tolerances are specified in `/app/SPEC.md`.
