"""Instrument design merit exploration from sampling, throughput, and RV information."""

from astropy.table import Table  # For reading FITS files
from astropy.io import fits  # For reading FITS files
import numpy as np  # For numerical operations
import matplotlib.pyplot as plt  # For plotting
#uis
from scipy.interpolate import InterpolatedUnivariateSpline as ius  # For spline interpolation
from scipy import constants  # For physical constants



# construct a wavelenght grid at exactly 1 km/s/pix
# Function to generate a logarithmic wavelength grid.
def get_magic_grid(wave0=360, wave1=920, dv_grid=0.5):
    """
    Generate a logarithmic wavelength grid.

    :param wave0: Starting wavelength.
    :param wave1: Ending wavelength.
    :param dv_grid: Velocity step in km/s.
    :return: Logarithmic wavelength grid.
    """
    # Calculate the number of grid points based on the velocity step.
    len_magic = int(np.ceil(np.log(wave1 / wave0) * np.array(constants.c / 1000) / dv_grid))
    # Generate the logarithmic wavelength grid.
    magic_grid = np.exp(np.arange(len_magic) / len_magic * np.log(wave1 / wave0)) * wave0
    return magic_grid

flux_file = '/Users/eartigau/vroomm_py/vroomm_simu/data/lte04000-5.00-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits'
wave_file = '/Users/eartigau/vroomm_py/vroomm_simu/data/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits'

transmission = Table.read('vroomm-frontend-transmission.csv')

wave_trans = transmission['wavelength']
transmission = transmission['transmission']

flux = fits.getdata(flux_file)
wave = fits.getdata(wave_file)/10.0

keep = (wave>350)*(wave<1000)
wave = wave[keep]
flux = flux[keep]
# express flux in photons, not erg
flux*=wave
flux/=np.nanmedian(flux)

wave_magic = get_magic_grid()
flux_magic = ius(wave, flux)(wave_magic)
flux_magic*=1e4

tbl = Table.read('/Users/eartigau/vroomm_py/data/VROOMM_orders_centroids.txt',format='ascii')
tbl['XPIX']*=(1000/12.0) # convert mm to pixels

waves = np.array([0.4,0.525,0.75])
steps = np.zeros(len(waves))

uord = np.unique([int(i) for i in tbl['ORDER'].data])

wave_ord = np.zeros(len(uord))
step_ord = np.zeros(len(uord))


for iord, ord in enumerate(uord):
    g = (tbl['ORDER'] == ord)
    step = np.mean((tbl[g]['WAVELENGTH']/np.gradient(tbl[g]['WAVELENGTH'],tbl[g]['XPIX'])))/constants.c*1000
    wave_ord[iord] = np.mean(tbl[g]['WAVELENGTH'])
    step_ord[iord] = step
    print(f'Order {ord} : step = {step:.3f} km/s/pix')


# Plot the flux
plt.figure(figsize=(10, 5))

dv = np.arange(-10,10,0.5)
# Gaussian profile
def gaussian_box(x, mu, fwhm,expo = 3.0):
    sigma = fwhm / (2 * (2 * np.log(2)**(1/expo)))**(1/expo)
    g = np.exp(-np.abs(x - mu) ** expo / (2 * sigma ** expo))
    g /= np.sum(g)
    return g

cuts = [0.4, 0.525, 0.75, 0.920]



g1 = (wave_magic>cuts[0])*(wave_magic<cuts[1])
g2 = (wave_magic>cuts[1])*(wave_magic<cuts[2])
g3 = (wave_magic>cuts[2])*(wave_magic<cuts[3])

trans_1 = np.mean(transmission[(wave_trans>cuts[0])*(wave_trans<cuts[1])])
trans_2 = np.mean(transmission[(wave_trans>cuts[1])*(wave_trans<cuts[2])])
trans_3 = np.mean(transmission[(wave_trans>cuts[2])*(wave_trans<cuts[3])])



fws = np.arange(2.5,5.5,0.5)

err = np.zeros((len(fws),3))

steps =  np.polyval(np.polyfit(wave_ord,step_ord,2),[0.4,0.55,0.75])

for i, fw in enumerate(fws):
    # Generate the Gaussian profile
    gg = gaussian_box(dv, 0, fw)

    noise = np.sqrt(flux_magic)

    flux_magic2 = np.convolve(flux_magic, gg, mode='same')

    grad = np.gradient(flux_magic2)/np.gradient(wave_magic)*wave_magic/constants.c

    err1 = np.sqrt(1/np.sum(grad[g1]**2/noise[g1]**2))
    err2 = np.sqrt(1/np.sum(grad[g2]**2/noise[g2]**2))
    err3 = np.sqrt(1/np.sum(grad[g3]**2/noise[g3]**2))
    print(f'{fw:.1f} km/s: {err1:.2f} m/s, {err2:.2f} m/s, {err3:.2f} m/s')
    err[i,0] = err1
    err[i,1] = err2
    err[i,2] = err3

fws_pix0 = fws/steps[0]
fws_pix1 = fws/steps[1]
fws_pix2 = fws/steps[2]
fit1 = np.polyfit(3.0/fws_pix0-1, 1/err[:,0]**2, 1)
fit2 = np.polyfit(3.0/fws_pix1-1, 1/err[:,1]**2, 1)
fit3 = np.polyfit(3.0/fws_pix2-1, 1/err[:,2]**2, 1)

print(f'fit1: {fit1[0]:.2f} {fit1[1]:.2f}')
print(f'fit2: {fit2[0]:.2f} {fit2[1]:.2f}')
print(f'fit3: {fit3[0]:.2f} {fit3[1]:.2f}')


plt.plot(fws_pix0,1/err[:,0]**2,'o-', label='350-450 nm')
plt.plot(fws_pix1,1/err[:,1]**2,'o-', label='450-600 nm')
plt.plot(fws_pix2,1/err[:,2]**2,'o-', label='600-920 nm')
plt.xlabel('3.0 km/s/FWHM')
plt.ylabel('Efficiency')
plt.title('Efficiency vs FWHM')
plt.legend()
plt.tight_layout()
plt.savefig('vroomm/efficiency.png', dpi=300)
plt.show()


print()

# L'efficacité va comme :

print('Domaine #1 : 350-450 nm')
print('Domaine #2 : 450-600 nm')
print('Domaine #3 : 600-920 nm')
print()
print(f'Mean pixel size for domain #1 : {steps[0]:.2f} km/s/pix')
print(f'Mean pixel size for domain #1 : {steps[1]:.2f} km/s/pix')
print(f'Mean pixel size for domain #1 : {steps[2]:.2f} km/s/pix')
print()

print(f'effic_1 = {fit1[0]:.2f} * (3.0/FWHM_1-1) + {fit1[1]:.2f}')
print(f'effic_2 = {fit2[0]:.2f} * (3.0/FWHM_2-1) + {fit2[1]:.2f}')
print(f'effic_3 = {fit3[0]:.2f} * (3.0/FWHM_3-1) + {fit3[1]:.2f}')
print()

print('FWHM expressed in pixels [12µm]')
print()

print('Fonction de mérite globale : ')
print(f'effic_1*throughput_1 + effic_2*throughput_2 + effic_3*throughput_3')