import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.signal import convolve
from scipy.interpolate import InterpolatedUnivariateSpline as ius
from scipy.constants  import c


"""
For a given resolution, what is the efficiency loss due to the fiber size?
This uses the LHS1140 template from the ESPRESSO pipeline and the NIRPS HE fiber size assuming
the use of a 2x slicer with a rectangular fiber.

This is compared with an octogonal fiber with the same area as the rectangular fiber.
"""

template = Table.read('/Users/eartigau/vroomm_simu/data/Template_LHS1140_tc_ESPRESSO.fits')

# expressed in arcsec
fib_size = [2.0,2.25,2.5,2.75,3.0,3.25,3.5,3.75]
fib_eff = [0.65,0.695,0.73,0.77,0.80,0.82,0.85,0.88]
loss_square = 0.85

fib_size = np.array(fib_size)
fib_eff = np.array(fib_eff)

f_telescope =8.0
f_fiber = 2.3
pix_scale = 12.5 # microns
scale_ratio = f_telescope/f_fiber
vsini = 10.0 # km/s

fib_size_kms = fib_size*62.5/scale_ratio/pix_scale # expressed in km/s at F/# of telescope
square_size_kms = fib_size*62.5/scale_ratio/pix_scale/2 # expressed in km/s at F/# of telescope

fiber_core_mu = fib_size*62.5/scale_ratio

wave1 = 650
wave2 = 700
template = template[template['wavelength']>wave1]
template = template[template['wavelength']<wave2]

# sampling = 250 m/s
flux = np.array(template['flux'])
flux[~np.isfinite(flux)] = 1.0

wave = np.array(template['wavelength'])

gg = (wave>669)*(wave<670)
q0 = np.sum(np.gradient(flux)[gg]**2)

qq1s = np.zeros(len(fib_size))
qq2s = np.zeros(len(fib_size))

resolution1 = (c/1000)/fib_size_kms*1.15
resolution2 = (c/1000)/square_size_kms*1.15

for ifib in range(len(fib_size)):
    w1 = int(fib_size_kms[ifib]*4.0)
    ker1 = np.sin(np.pi*(np.arange(w1)+0.5)/w1)
    ker1/=np.sum(ker1)

    w2 = int(square_size_kms[ifib]*4.0)

    ker2 =np.exp(-( (np.arange(w2+4)-(w2+4)/2.0+0.5)/(w2/2.0) )**20 )
    ker2/=np.sum(ker2)

    flux1 = convolve(flux, ker1, mode='same')
    flux2 = convolve(flux, ker2, mode='same')

    qq1s[ifib] = np.sum(np.gradient(flux1)[gg] ** 2)
    qq2s[ifib] = np.sum(np.gradient(flux2)[gg] ** 2)

qq1s*=fib_eff
qq2s*=fib_eff

q0 = qq2s[0].copy()
qq1s = qq1s/q0
qq2s = qq2s/q0

inirps = 1

fig, ax = plt.subplots(1,2, figsize = [8,4], sharey = True)
ax[0].plot(fib_size,qq1s,'g-', label = 'circular fiber')
ax[0].plot(fib_size,qq2s,'r-', label = 'rectangular fiber')
ax[0].plot(fib_size[inirps],qq2s[inirps],'ro', label = 'NIRPS HE')


ax[0].plot(fib_size,qq2s*loss_square,'c-', label = 'rectangular fiber ({:.2f} efficiency)'.format(loss_square))
ax[0].plot(fib_size[inirps],qq2s[inirps]*loss_square,'co', label = 'NIRPS HE')

ax[0].set_xlabel('Fiber size on sky (arcsec)')
ax[0].set_ylabel('Overall efficiency')
ax[0].legend(fontsize = 8)


ax[1].plot(resolution1,qq1s,'g-', label = 'circular fiber')
ax[1].plot(resolution2,qq2s,'r-', label = 'rectangular fiber')
ax[1].plot(resolution2,qq2s*loss_square,'c-', label = 'rectangular fiber ({:.2f} efficiency)'.format(loss_square))


ax[1].plot(resolution2[inirps],qq2s[inirps],'ro', label = 'NIRPS HE')
ax[1].plot(resolution2[inirps],qq2s[inirps]*loss_square,'co', label = 'NIRPS HE')


ax[1].set_xlabel('Resolution ($\lambda/\Delta\lambda$)')
ax[0].grid(color = 'grey', alpha = 0.5)
ax[1].grid(color = 'grey', alpha = 0.5)

#ax[1].legend()
ax[0].set_ylim([0.4,1.05])
plt.tight_layout()
plt.savefig('efficiency_loss.pdf', dpi = 300)
plt.show()


