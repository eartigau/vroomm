from astroquery.gaia import Gaia
from astropy.table import Table
import os
from astroquery.simbad import Simbad
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt

gaia_table_name = "/Users/eartigau/vroomm_simu/big_gaia_table.csv"

# on sauvegarde la table localement pour s'économiser le download la prochaine
# fois
if not os.path.isfile(gaia_table_name):
    QUERY = "SELECT * FROM gaiadr3.gaia_source WHERE parallax > 7 AND parallax < 7.5 AND ra > 240 AND ra < 255 AND " \
            "dec > -30 AND dec < -20 "
    job = Gaia.launch_job_async(QUERY)

    gaia_table_big = job.get_results()

    gaia_table_big.write(gaia_table_name,overwrite = True)

gaia_table_big = Table.read(gaia_table_name)


dist = 1000/gaia_table_big['parallax']*3.0857e13 # distance in km

pmra_err = gaia_table_big['pmra_error']

velo_err = pmra_err/(206254*1000)*dist/(86400*365.24) # in km/s

rp = gaia_table_big['phot_rp_mean_mag']

# 10 min, R=14, 50m/s

mag_simu = np.arange(10,23)
rv_simu =  np.sqrt(10**( (mag_simu-14)/2.5))*0.05

plt.plot(rp,velo_err*1000,'k.', alpha = 0.15)
plt.plot(mag_simu,rv_simu*1000,'r-', alpha = 0.8, label = 'VROOMM - 10 min', linewidth = 2)
plt.plot(mag_simu,rv_simu/np.sqrt(6)*1000,'-', color = 'orange', alpha = 0.8, label = 'VROOMM - 1 h', linewidth = 2)
plt.xlabel('R magnitude')
plt.ylabel('GAIA sky-planet velocity error [m/s]')
plt.xlim([10,20])
plt.ylim([10,3000])
# log on y axis
plt.yscale('log')
plt.title('Upper Scorpius line of sight & distance')
plt.legend()
plt.tight_layout()
plt.savefig('vroomm/plx_rv_vroomm.png',dpi = 300)
plt.show()