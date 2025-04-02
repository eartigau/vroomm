import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.io import fits
import numpy as np
from etienne_tools import lowpassfilter, doppler
# import ius
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.optimize import curve_fit
from tqdm import tqdm
from scipy.signal import convolve
from astropy.constants import c
# get rotbroad
from PyAstronomy.pyasl import rotBroad

def gaussian(x, mu, sig, amp):

    return 1-amp*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))*np.sign(sig)

exptime = 1.0 #s
sampling = 1.0 #km/s/pix

r_mag = 7

wave1 = 600
wave2 = 850

r_primary = 160  # cm
wave0 = 6204.29  # A
resolution = 120000  # lambda/delta_lambda
npix = c.value/resolution/1000 # pixel per resolution element
effic = 0.1  # 10% efficiency

dark_current = 1.6e-4  # e-/s/pixel
spatial_extent = 10 # pixels in height
dark_current *= exptime
dark_current *= spatial_extent

area = np.pi * (r_primary / 2) ** 2
zp = 2.49767e-9  # (erg/cm2/s/A)
zp = zp * area  # (erg/s/A)

res_element = wave0 / resolution

zp = zp * res_element * exptime * 10 ** (-r_mag / 2.5) * effic  # erg

energ = 1.602176634e-12 * (12398.418 / wave0)

nphot = zp / energ  # photons

nphot_per_spec_pix = nphot/npix

nphot_per_pix = nphot/spatial_extent

print('Number of photons per pixel: {:.1f}, mag = {}'.format(nphot_per_pix,r_mag))

