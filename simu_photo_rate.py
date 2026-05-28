"""Quick photon-rate and integrated SNR scaling calculator for VROOMM scenarios."""

# Import necessary libraries
import matplotlib.pyplot as plt  # For plotting
from astropy.table import Table  # For handling tabular data
from astropy.io import fits  # For reading FITS files
import numpy as np  # For numerical operations
def lowpassfilter(input_vect, width=101):
    """Low-pass filter via running NaN-median spline, handling NaN values."""
    index = np.arange(len(input_vect))
    xmed = []
    ymed = []
    for i in np.arange(-width // 2, len(input_vect) + width // 2, width // 4):
        low_bound = max(i, 0)
        high_bound = min(i + int(width), len(input_vect) - 1)
        pixval = index[low_bound:high_bound]
        if len(pixval) < 3:
            continue
        if np.max(np.isfinite(input_vect[pixval])) == 0:
            continue
        xmed.append(np.nanmean(pixval))
        ymed.append(np.nanmedian(input_vect[pixval]))
    xmed = np.array(xmed, dtype=float)
    ymed = np.array(ymed, dtype=float)
    if len(xmed) < 3:
        return np.zeros_like(input_vect) + np.nan
    if len(xmed) != len(np.unique(xmed)):
        xmed2 = np.unique(xmed)
        ymed2 = np.zeros_like(xmed2)
        for i in range(len(xmed2)):
            ymed2[i] = np.mean(ymed[xmed == xmed2[i]])
        xmed = xmed2
        ymed = ymed2
    spline = ius(xmed, ymed, k=2, ext=3)
    return spline(np.arange(len(input_vect)))


def doppler(wave, v):
    """Relativistic Doppler shift. v in m/s."""
    v = np.array(v)
    wave = np.array(wave)
    return wave * np.sqrt((1 - v / c.value) / (1 + v / c.value))
from scipy.interpolate import InterpolatedUnivariateSpline as ius  # For spline interpolation
from scipy.optimize import curve_fit  # For curve fitting
from tqdm import tqdm  # For progress bars
from scipy.signal import convolve  # For signal convolution
from astropy.constants import c  # Speed of light constant
from PyAstronomy.pyasl import rotBroad  # For rotational broadening

# Define a Gaussian function for curve fitting
def gaussian(x, mu, sig, amp):
    """
    Returns a Gaussian function with the given parameters.

    Parameters:
    x (array): The input x values.
    mu (float): The mean of the Gaussian.
    sig (float): The standard deviation of the Gaussian.
    amp (float): The amplitude of the Gaussian.

    Returns:
    array: The Gaussian function evaluated at x.
    """
    return 1 - amp * np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) * np.sign(sig)

# Simulation parameters
exptime = 0.25  # Exposure time in seconds
t_tot = 600
sampling = 1.0  # Sampling in km/s/pixel

# Apparent magnitude in the r-band
r_mag = 10

# Wavelength range for the simulation (in nm)
wave1 = 600  # Start wavelength
wave2 = 850  # End wavelength

# Telescope and instrument parameters
r_primary = 160  # Telescope primary mirror radius in cm
wave0 = 6204.29  # Central wavelength in Angstroms
resolution = 80000  # Spectral resolution (lambda/delta_lambda)
npix = c.value / resolution / 1000  # Pixels per resolution element (in km/s/pixel)
effic = 0.15  # Efficiency (10%)

# Detector parameters
dark_current = 1.6e-4  # Dark current in e-/s/pixel
spatial_extent = 10  # Spatial extent in pixels (height of the spectrum)
spectral_extent = 3
dark_current *= exptime  # Scale dark current by exposure time
dark_current *= spatial_extent  # Scale dark current by spatial extent

# Calculate the telescope's collecting area
area = np.pi * (r_primary / 2) ** 2  # Mirror area in cm^2

# Zero-point flux in erg/cm^2/s/A
zp = 2.49767e-9  # Zero-point flux for a magnitude 0 star
zp = zp * area  # Scale by the telescope's collecting area (erg/s/A)

# Resolution element in Angstroms
res_element = wave0 / resolution

# Scale zero-point flux by resolution, exposure time, magnitude, and efficiency
zp = zp * res_element * exptime * 10 ** (-r_mag / 2.5) * effic  # erg

# Calculate the energy of a photon at the central wavelength
energ = 1.602176634e-12 * (12398.418 / wave0)  # Energy in erg

# Calculate the total number of photons collected
nphot = zp / energ  # Total number of photons

# Calculate the number of photons per spectral pixel
nphot_per_spec_pix = nphot / npix

# Calculate the number of photons per spatial pixel
nphot_per_pix = nphot / spatial_extent

# Print the result
print('Number of photons per pixel: {:.2f}, mag = {}'.format(nphot_per_pix, r_mag))
snr = np.sqrt(nphot_per_pix/exptime*t_tot*spatial_extent*3*spectral_extent)  # SNR calculation
print('SNR: {:.2f} in {:.2f} min'.format(snr, t_tot/60))
