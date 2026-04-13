# Script Documentation

This file is intended as a practical reference for graduate students.

## etc_vroomm.py

Purpose:
- Compute SNR per resolution element vs r magnitude for one or more exposure times.
- Compare baseline OMM, EMCCD-like readout behavior, and a larger telescope scaling case.

How to run:
- python etc_vroomm.py

Main parameters to edit first:
- exptime in the loop near the top
- resolution
- effic
- ron
- dark_current

Outputs:
- vroomm_snr.png

Quick interpretation:
- Horizontal lines mark rough RV precision regimes (1, 5, 100 m/s).
- Curves crossing these lines identify target magnitude limits for each setup.

## simu_photo_rate.py

Purpose:
- Back-of-the-envelope photon and SNR scaling for a fixed setup.

How to run:
- python simu_photo_rate.py

Main parameters to edit first:
- r_mag
- exptime
- t_tot
- resolution
- effic

Outputs:
- Console summary:
	- photons per pixel
	- integrated SNR over total time

Quick interpretation:
- Use this script first when testing sensitivity to magnitude and throughput assumptions.

## simu_faint.py

Purpose:
- Simulate noisy spectra for faint targets and retrieve RV through CCF fitting.

How to run:
- python simu_faint.py

Required inputs:
- template FITS spectrum (currently hard-coded absolute path)
- CCF mask FITS (currently hard-coded absolute path)
- local etienne_tools module (lowpassfilter, doppler)

Main parameters to edit first:
- r_mag
- exptime
- wave1, wave2
- nmc
- vsini

Outputs:
- ccf_rmag_<mag>_exptime_<seconds>.pdf

Quick interpretation:
- Inspect RV scatter and CCF width from the Monte Carlo ensemble.
- For fair comparisons, vary one physical parameter at a time.

## effic_loss.py

Purpose:
- Quantify resolution/efficiency trade-offs for fiber geometry choices.

How to run:
- python effic_loss.py

Required inputs:
- Template_LHS1140_tc_ESPRESSO.fits (absolute path in current code)

Main parameters to edit first:
- fib_size
- fib_eff
- loss_square
- wave1, wave2

Outputs:
- efficiency_loss.pdf

Quick interpretation:
- Left panel: efficiency trend vs on-sky fiber size.
- Right panel: same trend mapped to spectral resolution.

## optimisation_design.py

Purpose:
- Build a simple merit trend linking line-spread function width, wavelength domain, and RV information content.

How to run:
- python optimisation_design.py

Required inputs:
- PHOENIX flux FITS
- PHOENIX wavelength FITS
- vroomm-frontend-transmission.csv
- VROOMM_orders_centroids.txt

Main parameters to edit first:
- cuts (wavelength domains)
- fws (trial FWHM values)
- get_magic_grid velocity step

Outputs:
- Currently writes vroomm/efficiency.png
- Prints fitted linearized efficiency laws per domain

Quick interpretation:
- Resulting slopes and intercepts are useful for first-pass design scoring.

## plx_rv_vroomm.py

Purpose:
- Compare Gaia tangential velocity precision with empirical VROOMM RV scaling.

How to run:
- python plx_rv_vroomm.py

Inputs:
- Local Gaia table cache, or query Gaia archive if cache does not exist

Main parameters to edit first:
- Gaia ADQL query box
- rv_simu normalization relation

Outputs:
- Currently writes vroomm/plx_rv_vroomm.png

Quick interpretation:
- Helps identify where RV follow-up precision can match or exceed sky-plane constraints.

## scripts/generate_public_plots.py

Purpose:
- Generate two public-safe figures without private FITS inputs.

How to run:
- python scripts/generate_public_plots.py

Outputs:
- figures/vroomm_snr_public.png
- figures/vroomm_photon_snr_grid.png

When to use:
- First test after cloning.
- Sanity check that your environment is configured correctly.
