import numpy as np
import matplotlib.pyplot as plt

r_mag = np.arange(5,28)

frac_bandpass_guinding = 4
guiding_throughput = 0.04


i=0
colors = ['red','orange','green','blue']
fig, ax = plt.subplots(layout='constrained')
for exptime in [1800]:
    #exptime = 300.0 # s

    r_primary = 160 # cm
    wave0 = 6204.29 # A
    resolution = 120000 # lambda/delta_lambda
    effic = 0.1 # 10% efficiency
    spatial_exten = 2*10# pixels
    ron = 3.0 # readout noise
    dark_current_4m = 2.9/3600 # e-/s/pixel
    dark_current = 1.6e-4 # e-/s/pixel

    area = np.pi*(r_primary/2)**2
    zp = 2.49767e-9#	(erg/cm2/s/A)
    zp = zp*area # (erg/s/A)

    res_element = wave0/resolution

    zp = zp*res_element*exptime*10**(-r_mag/2.5)*effic # erg

    energ = 1.602176634e-12*(12398.418/wave0)

    nphot = zp/energ # photons

    print(nphot)

    if exptime == 60:
        nphot_per_sec = nphot/exptime

        guiding_photons = nphot_per_sec*(resolution/frac_bandpass_guinding)*guiding_throughput/effic
        for imag in range(len(r_mag)):
            print('guiding photons per second {:.0f}, mag = {:.0f}'.format(guiding_photons[imag],r_mag[imag]))

    snr = nphot/np.sqrt(nphot+ron**2*spatial_exten+dark_current*spatial_exten*exptime)
    plt.plot(r_mag,snr,label = str(exptime)+' s - OMM', color = colors[i])

    snr = nphot/np.sqrt(nphot+dark_current*spatial_exten*exptime)
    plt.plot(r_mag,snr,label = str(exptime)+' s - OMM, EMCCD', color = colors[i],linestyle = '--')

    snr2 = nphot*6/np.sqrt(nphot*6+ron**2*spatial_exten+dark_current_4m*spatial_exten*exptime)
    plt.plot(r_mag,snr2,label = str(exptime)+' s - 4-m telescope', color = colors[i],linestyle = 'dashdot')


    i+=1
ax.set_xlabel('r magnitude')
ax.set_ylabel('SNR per resolution element')
#ax.grid(color = 'grey',alpha = 0.5)
# get log on y axis
ax.set_yscale('log')
ax.plot([5,30],[1,1],'k:',alpha = 0.5, label = '100 m/s')
ax.plot([5,30],[200,200],'k--',alpha = 0.5, label = '1 m/s')
ax.plot([5,30],[40,40],'k-.',alpha = 0.5, label = '5 m/s')

ax.set(ylim = [0.1,600])
ax.set(xlim = [5,26])
plt.legend()
ax.set(title = 'VROOMM SNR')
plt.savefig('vroomm_snr.png')
plt.show()
