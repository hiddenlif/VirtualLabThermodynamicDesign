# Virtual Lab: thermodynamic design

Choose calorimetry experiments and use their simulated measurements to certify the operating response of two synthetic protein lots. Experiment choices affect both the observations and irreversible stock consumption.

**Task type:** Interactive simulated experimentation. **Source version:** 2.0.2.

All scientific specimens, measurements and instrument conditions are synthetic. These examples are not validated biological assays or clinical decision tools.

## Contents

- `TASK.md`: scientific task and required deliverable.
- `app/`: public scientific model, documentation and disclosed practice inputs where supplied.
- `Dockerfile`: the task's original dependency environment.

## Use

Build the environment with `docker build -t virtual-lab-thermodynamics .`.

Implement `app/controller.py` using the JSON-lines contract in `TASK.md` and `app/SPEC.md`. `app/campaign.py` exposes `Campaign.start()` and `Campaign.batch(experiments)` for the disclosed practice case. Its command-line entry point only loads and initializes a campaign; the complete controller evaluation runs in the private verifier. The public simulator returns action-dependent measurements and stock state.

## Execution evidence

Three historical reference runs passed, each accepting all three campaigns. A no-op scored zero. These are local reference results, not blind solver results or platform approval.

The historical evidence concerns the source version. This publication separates public inputs from private evaluation material and cleans documentation. No new task execution or model trial was performed for this publication, and the cleaned repository has not been revalidated end to end.

## Reference and verifier

Reference implementations, held-out cases, expected results and grading code are kept in the private repository `science-task-solutions`, under `VirtualLabThermodynamicDesign/`. They are excluded from this public repository. The public practice inputs and worked examples remain public by design.
