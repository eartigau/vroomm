# Graduate Student Onboarding Guide

This guide is written for new members joining the VROOMM simulation effort.

## Week 1 Checklist

1. Create environment and run one script end-to-end.
2. Produce two reproducible figures from scripts/generate_public_plots.py.
3. Read docs/SCRIPTS.md and identify one script you will own.
4. Write down all hard-coded assumptions you changed.

## Recommended Learning Path

1. Start with etc_vroomm.py to build intuition for SNR scaling.
2. Use simu_photo_rate.py for quick observing-time intuition.
3. Move to simu_faint.py once you are ready for Monte Carlo RV retrieval.
4. Use effic_loss.py and optimisation_design.py for instrument trade studies.
5. Use plx_rv_vroomm.py for science-context comparisons with Gaia.

## Naming and Output Conventions

- Save plots with descriptive names containing key parameters.
- Keep outputs in figures/ whenever practical.
- If a script writes to a non-existing subfolder, create the folder first or update the output path.

## Reproducibility Habits

- Record exact script name and edited constants in your note.
- Record commit hash with each figure used in a report.
- Change one physical assumption at a time when exploring sensitivity.

## Scientific Caution

- Many scripts are prototype-level and optimized for speed of exploration.
- Units are mixed across files; verify unit consistency before interpretation.
- Absolute paths indicate local development history, not portable workflow.

## First Refactor Tasks (Good Starter Contributions)

1. Replace one absolute input path by a path relative to the repository.
2. Add argparse to one script so key parameters can be passed at runtime.
3. Add one validation plot that catches common configuration mistakes.
4. Move one repeated constant block into a shared helper function.

## Asking For Help

When asking for review, include:
- objective of the run
- exact command executed
- edited constants
- one plot and one sentence explaining the result
