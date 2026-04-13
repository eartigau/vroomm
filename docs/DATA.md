# Data and Path Policy

Several scripts use machine-specific absolute paths (for example under `/Users/...`).

## Why this matters

Absolute paths break reproducibility on GitHub and CI.

## Migration target

Use a local repository structure:

- `data/templates/`
- `data/masks/`
- `data/models/`
- `data/catalogs/`
- `data/design/`

## Minimal migration pattern

Replace:

```python
Table.read('/Users/name/project/data/file.fits')
```

With:

```python
from pathlib import Path
DATA = Path(__file__).resolve().parent / 'data'
Table.read(DATA / 'templates' / 'file.fits')
```

## Large files

If files are large:
- keep lightweight examples in git
- host heavy products externally (institution server or object storage)
- document download instructions in this file

## Optional enhancement

Add a single environment variable for data location:

- `VROOMM_DATA_DIR`

Then resolve paths with fallback to `./data`.

## Script-Specific Inputs

This section maps each script to the external files it expects.

- simu_faint.py
	- Template spectrum FITS
	- CCF mask FITS
	- local etienne_tools module

- effic_loss.py
	- Template_LHS1140_tc_ESPRESSO.fits

- optimisation_design.py
	- PHOENIX flux FITS
	- PHOENIX wavelength FITS
	- vroomm-frontend-transmission.csv
	- VROOMM_orders_centroids.txt

- plx_rv_vroomm.py
	- Gaia cache table path (or online query fallback)

## Suggested Portable Layout

If you are setting up from scratch, this layout works well:

- data/templates/
- data/masks/
- data/models/
- data/catalogs/
- data/design/

## Common Errors and Fixes

- FileNotFoundError on a FITS input
	- Confirm the file is present under your local data folder.
	- Confirm the script path points to that folder.

- Script writes into a missing folder like vroomm/
	- Create the folder or update output paths to figures/.

- Astroquery timeout or archive unavailable
	- Use a local cached table for development runs.

- ImportError for etienne_tools
	- Add the local tools module to your Python path or vendor required functions.
