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

exptime = 3600#s
sampling = 1.0 #km/s/pix

r_mag = 18

wave1 = 600
wave2 = 850

vsini = 0.1 # km/s
wave_vsini = (1+np.arange(np.round(vsini*2+3))/(c.value/1000))*wave1
sp_vsini = np.zeros(len(wave_vsini))
sp_vsini[len(sp_vsini)//2] = 1.0
ker = rotBroad(wave_vsini, sp_vsini, 0.6, float(vsini),edgeHandling="firstlast")

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

npot_per_pix = nphot/npix

tbl = Table.read('/Users/eartigau/vroomm_simu/data/Template_LHS1140_tc_ESPRESSO.fits')
tbl = tbl[:4*int(len(tbl)//4)]
sampling = int(c.value/np.nanmedian(tbl['wavelength']/np.gradient(tbl['wavelength'])))/1000. # pix/km/s

mask = Table.read('/Users/eartigau/vroomm_simu/data/GL846_tc_full.fits')

mask = mask[(mask['ll_mask_s']>wave1) & (mask['ll_mask_s']<wave2)]
#mask = mask[mask['w_mask']>0]
wave_mask = np.array(mask['ll_mask_s'])
weight_mask = np.array(mask['w_mask'])


w = np.array(tbl['wavelength'])
f = np.array(tbl['flux'], dtype =float)
f[~np.isfinite(f)] = 1
f = (f[::4]+f[1::4]+f[2::4]+f[3::4])/4 # we get to 1 km/s/pix
w = (w[::4]+w[1::4]+w[2::4]+w[2::4])/4

f = np.convolve(f,ker,mode = 'same')
f/=lowpassfilter(f,1000)

g = (w>wave1) & (w<wave2)
w = w[g]
f = f[g]

nmc = 30
rv = np.zeros(nmc)
fw = np.zeros(nmc)

offset = np.arange(-150, 150, 1)
ccfs = np.zeros((nmc,len(offset)))
for ii in tqdm(range(nmc)):
    signal = np.random.poisson(f*npot_per_pix+dark_current)-dark_current

    spl = ius(w,signal, k=1, ext=1)

    ccf = np.zeros(len(offset))
    for i in range(len(offset)):
        mask2 = doppler(wave_mask,offset[i]*1000)
        ccf[i] = np.nansum(spl(mask2)*weight_mask)
    ccf/=np.nanmedian(ccf)
    ccfs[ii] = ccf

    # mu, sig, amp
    ccf2 = convolve(ccf, np.ones(5)/5, mode='same')
    ccf2[0:5] = 1.0
    ccf2[-5:] = 1.0
    p0 = offset[np.argmin(ccf2)], 3, 1-ccf2[np.argmin(ccf2)]
    fit, _ = curve_fit(gaussian,offset,ccf,p0=p0)
    rv[ii] = fit[0]*2.45
    fw[ii] = fit[1]*2.45


n1, med, p1 = np.percentile(ccfs, [16, 50, 84],axis=0)

fig, ax = plt.subplots(nrows = 2, ncols = 1, sharex = False, figsize = (16 ,8))
ax[0].plot(w,signal,alpha = 0.5, label = 'Signal + noise')
ax[0].plot(w,f*npot_per_pix,alpha = 0.9, color = 'k', label = 'Input signal')
ax[0].set(xlim=[700,701])
ax[0].set(xlabel = 'Wavelength (nm)')
ax[0].set(ylabel = 'Flux (photons/pix)')
ax[1].set(xlabel = 'RV (km/s)', ylabel = 'CCF')
ax[0].set(ylim = [-npot_per_pix,npot_per_pix*3])

ax[1].plot(offset,ccf, label = 'sample CCF',color = 'green')
ax[1].fill_between(offset,n1,p1, alpha = 0.2, color = 'green',  label = '1-$\sigma$ envelope')

snr = npot_per_pix/(np.nanstd(signal  - np.roll(signal,1))/np.sqrt(2))

ax[1].plot(offset,gaussian(offset,*fit),color = 'red', alpha = 0.5, label = 'Gaussian fit')
ax[1].text(0.6,0.05,'RV = {:.3f}$\pm${:.3f} km/s\nfwhm = {:.2f}$\pm${:.2f} km/s'.format(np.nanmean(rv), np.nanstd(rv),
                                                                                    np.nanmean(fw), np.nanstd(fw)),
           transform=ax[1].transAxes,  bbox = dict(boxstyle='round', facecolor='wheat', alpha=0.5))
ax[1].legend()
ax[0].legend()
ax[1].grid(color = 'k', alpha = 0.2)
ax[0].set(title = 'Rmag = {:.1f}, exptime = {}s, SNR = {:.2f}'.format(r_mag, exptime,snr))
ax[1].set(xlim = [np.min(offset),np.max(offset)])
plt.tight_layout()
outname = 'ccf_rmag_{:.1f}_exptime_{:.0f}.pdf'.format(r_mag, exptime)
plt.savefig(outname)
plt.show()

