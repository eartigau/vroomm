# VROOMM Simulation Toolkit

Design and performance simulation scripts for VROOMM, a high-resolution spectrograph concept for the Mont Megantic Observatory.

This repository is aimed at exploratory instrument studies and is now documented for new graduate students joining the project.

## What You Can Do With This Repo

- Estimate SNR as a function of magnitude and exposure time.
- Build first-order photon budget and detectability intuition.
- Explore RV precision trends from CCF simulations.
- Study efficiency trade-offs driven by fiber geometry and sampling.
- Compare expected VROOMM RV performance against Gaia astrometric constraints.

## Fast Start For New Students

### 1) Create and activate environment (conda)

```bash
conda create -n vroomm python=3.12 -y
conda activate vroomm
pip install -r requirements.txt
```

### 2) Run a fully reproducible script first

```bash
python scripts/generate_public_plots.py
```

Expected outputs:
- figures/vroomm_snr_public.png
- figures/vroomm_photon_snr_grid.png

### 3) Read project docs in this order

1. docs/GRAD_STUDENT_GUIDE.md
2. docs/SCRIPTS.md
3. docs/DATA.md

## Code Map

- etc_vroomm.py
	- Analytic SNR curves versus r magnitude.
	- Output: vroomm_snr.png

- simu_photo_rate.py
	- Quick photon-rate and integrated SNR estimate.
	- Output: console-only summary values.

- simu_faint.py
	- Monte Carlo CCF simulation for faint targets.
	- Output: ccf_rmag_<mag>_exptime_<seconds>.pdf
	- Requires local template and mask FITS files.

- effic_loss.py
	- Fiber geometry versus efficiency/resolution trade study.
	- Output: efficiency_loss.pdf
	- Requires local template FITS file.

- optimisation_design.py
	- Merit-function exploration vs sampling and throughput.
	- Output: currently writes vroomm/efficiency.png
	- Requires local PHOENIX model and design tables.

- plx_rv_vroomm.py
	- Gaia proper-motion error compared to simple RV scaling law.
	- Output: currently writes vroomm/plx_rv_vroomm.png

- scripts/generate_public_plots.py
	- Public-safe, no-private-data plotting script.
	- Outputs figures in figures/

## License

MIT License. See LICENSE.
