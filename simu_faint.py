# Import necessary libraries
import matplotlib.pyplot as plt  # For plotting
from astropy.table import Table  # For handling tabular data
from astropy.io import fits  # For reading FITS files
import numpy as np  # For numerical operations
from etienne_tools import lowpassfilter, doppler  # Custom tools for filtering and Doppler shift
from scipy.interpolate import InterpolatedUnivariateSpline as ius  # For spline interpolation
from scipy.optimize import curve_fit  # For curve fitting
from tqdm import tqdm  # For progress bars
from scipy.signal import convolve  # For signal convolution
from astropy.constants import c  # Speed of light constant
from PyAstronomy.pyasl import rotBroad  # For rotational broadening

# Define a Gaussian function for curve fitting
def gaussian(x, mu, sig, amp):
    # Returns a Gaussian function with given mean (mu), standard deviation (sig), and amplitude (amp)
    return 1 - amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) * np.sign(sig)

# Simulation parameters
exptime = 3600  # Exposure time in seconds
sampling = 1.0  # Sampling in km/s/pixel
r_mag = 18.5  # Apparent magnitude in the r-band
wave1 = 600  # Start wavelength in nm
wave2 = 850  # End wavelength in nm
vsini = 0.1  # Rotational velocity in km/s

# Generate a wavelength grid for rotational broadening
wave_vsini = (1 + np.arange(np.round(vsini * 2 + 3)) / (c.value / 1000)) * wave1
sp_vsini = np.zeros(len(wave_vsini))  # Initialize spectrum
sp_vsini[len(sp_vsini) // 2] = 1.0  # Set central value to 1
ker = rotBroad(wave_vsini, sp_vsini, 0.6, float(vsini), edgeHandling="firstlast")  # Apply rotational broadening

# Telescope and instrument parameters
r_primary = 160  # Telescope primary mirror radius in cm
wave0 = 6204.29  # Central wavelength in Angstroms
resolution = 120000  # Spectral resolution (lambda/delta_lambda)
npix = c.value / resolution / 1000  # Pixels per resolution element
effic = 0.1  # Efficiency (10%)

# Detector parameters
dark_current = 1.6e-4  # Dark current in e-/s/pixel
spatial_extent = 10  # Spatial extent in pixels
dark_current *= exptime  # Scale dark current by exposure time
dark_current *= spatial_extent  # Scale by spatial extent

# Calculate zero-point flux and photon count
area = np.pi * (r_primary / 2) ** 2  # Mirror area in cm^2
zp = 2.49767e-9  # Zero-point flux in erg/cm^2/s/A
zp = zp * area  # Scale by mirror area
res_element = wave0 / resolution  # Resolution element in Angstroms
zp = zp * res_element * exptime * 10 ** (-r_mag / 2.5) * effic  # Scale by exposure time, magnitude, and efficiency
energ = 1.602176634e-12 * (12398.418 / wave0)  # Photon energy in erg
nphot = zp / energ  # Total number of photons
npot_per_pix = nphot / npix  # Photons per pixel

# Load template spectrum
tbl = Table.read('/Users/eartigau/vroomm_simu/data/Template_LHS1140_tc_ESPRESSO.fits')
tbl = tbl[:4 * int(len(tbl) // 4)]  # Trim table to a multiple of 4
sampling = int(c.value / np.nanmedian(tbl['wavelength'] / np.gradient(tbl['wavelength']))) / 1000.  # Sampling in km/s/pix

# Load mask for cross-correlation
mask = Table.read('/Users/eartigau/vroomm_simu/data/GL846_tc_full.fits')
mask = mask[(mask['ll_mask_s'] > wave1) & (mask['ll_mask_s'] < wave2)]  # Filter mask by wavelength range
wave_mask = np.array(mask['ll_mask_s'])  # Mask wavelengths
weight_mask = np.array(mask['w_mask'])  # Mask weights

# Process template spectrum
w = np.array(tbl['wavelength'])  # Wavelength array
f = np.array(tbl['flux'], dtype=float)  # Flux array
f[~np.isfinite(f)] = 1  # Replace non-finite values with 1
f = (f[::4] + f[1::4] + f[2::4] + f[3::4]) / 4  # Downsample to 1 km/s/pix
w = (w[::4] + w[1::4] + w[2::4] + w[2::4]) / 4  # Downsample wavelength array

# Apply rotational broadening and normalize
f = np.convolve(f, ker, mode='same')  # Convolve with rotational kernel
f /= lowpassfilter(f, 1000)  # Normalize using a low-pass filter

# Select wavelength range
g = (w > wave1) & (w < wave2)  # Boolean mask for wavelength range
w = w[g]  # Filter wavelength array
f = f[g]  # Filter flux array

# Monte Carlo simulation for radial velocity measurement
nmc = 30  # Number of Monte Carlo iterations
rv = np.zeros(nmc)  # Array to store radial velocities
fw = np.zeros(nmc)  # Array to store full-width at half-maximum (FWHM)
offset = np.arange(-150, 150, 1)  # RV offset range in km/s
ccfs = np.zeros((nmc, len(offset)))  # Array to store cross-correlation functions (CCFs)

for ii in tqdm(range(nmc)):  # Loop over Monte Carlo iterations
    # Simulate noisy signal
    signal = np.random.poisson(f * npot_per_pix + dark_current) - dark_current

    # Interpolate signal
    spl = ius(w, signal, k=1, ext=1)

    # Compute cross-correlation function (CCF)
    ccf = np.zeros(len(offset))
    for i in range(len(offset)):
        mask2 = doppler(wave_mask, offset[i] * 1000)  # Doppler-shift mask
        ccf[i] = np.nansum(spl(mask2) * weight_mask)  # Cross-correlation
    ccf /= np.nanmedian(ccf)  # Normalize CCF
    ccfs[ii] = ccf  # Store CCF

    # Fit Gaussian to CCF
    ccf2 = convolve(ccf, np.ones(5) / 5, mode='same')  # Smooth CCF
    ccf2[0:5] = 1.0  # Avoid edge effects
    ccf2[-5:] = 1.0
    p0 = offset[np.argmin(ccf2)], 3, 1 - ccf2[np.argmin(ccf2)]  # Initial guess for Gaussian fit
    fit, _ = curve_fit(gaussian, offset, ccf, p0=p0)  # Fit Gaussian
    rv[ii] = fit[0] * 2.45  # Convert to km/s
    fw[ii] = fit[1] * 2.45  # Convert FWHM to km/s

# Compute percentiles for CCF envelope
n1, med, p1 = np.percentile(ccfs, [16, 50, 84], axis=0)

# Plot results
fig, ax = plt.subplots(nrows=2, ncols=1, sharex=False, figsize=(16, 8))

# Plot signal and input spectrum
ax[0].plot(w, signal, alpha=0.5, label='Signal + noise')
ax[0].plot(w, f * npot_per_pix, alpha=0.9, color='k', label='Input signal')
ax[0].set(xlim=[700, 701], xlabel='Wavelength (nm)', ylabel='Flux (photons/pix)')
ax[0].set(ylim=[-npot_per_pix, npot_per_pix * 3])

# Plot CCF and Gaussian fit
ax[1].plot(offset, ccf, label='Sample CCF', color='green')
ax[1].fill_between(offset, n1, p1, alpha=0.2, color='green', label='1-$\sigma$ envelope')
snr = npot_per_pix / (np.nanstd(signal - np.roll(signal, 1)) / np.sqrt(2))  # Compute SNR
ax[1].plot(offset, gaussian(offset, *fit), color='red', alpha=0.5, label='Gaussian fit')
ax[1].text(0.6, 0.05, 'RV = {:.3f}$\pm${:.3f} km/s\nFWHM = {:.2f}$\pm${:.2f} km/s'.format(
    np.nanmean(rv), np.nanstd(rv), np.nanmean(fw), np.nanstd(fw)),
           transform=ax[1].transAxes, bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax[1].legend()
ax[0].legend()
ax[1].grid(color='k', alpha=0.2)
ax[0].set(title='Rmag = {:.1f}, exptime = {}s, SNR = {:.2f}'.format(r_mag, exptime, snr))
ax[1].set(xlim=[np.min(offset), np.max(offset)], xlabel='RV (km/s)', ylabel='CCF')

# Save and show plot
plt.tight_layout()
outname = 'ccf_rmag_{:.1f}_exptime_{:.0f}.pdf'.format(r_mag, exptime)
plt.savefig(outname)
plt.show()