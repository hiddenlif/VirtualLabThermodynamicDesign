# Virtual Lab Task: Calorimetry Experiment Design

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

On 2026-10-24 this published repository was run end to end: with the private reference installed in `app/`, the private grader accepted all three campaigns (reward 1, about 10 minutes); the unchanged starter scored 0. Three historical reference runs of the source version also passed. These are local reference results in Docker, not blind solver results or platform approval.

## Reference and verifier

Reference implementations, held-out cases, expected results and grading code are kept in the private repository `science-task-solutions`, under `VirtualLabThermodynamicDesign/`. They are excluded from this public repository. The public practice inputs and worked examples remain public by design.
