"""plot the [6.4] - [9.3] figure in the paper"""


import os
import math
import numpy as np
from scipy.interpolate import splrep, splev
import ss_readfits
import matplotlib.pyplot as plt
from scipy.integrate import simps
from astropy.stats import median_absolute_deviation
import matplotlib
from matplotlib import rc

matplotlib.rcParams.update({'font.size': 35})
rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)


g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")

f, axarr = plt.subplots(3, sharex=True, figsize = (20, 45))
plt.xlabel('[6.4] - [9.3]')



c6 = np.linspace(6.25,6.55,11)
c9 = np.linspace(9.1, 9.5, 8)
base = 100**0.2
s6 = 96.5
s9 = 45.7

ms = 1.98892e30 #kg
G = 6.67e-11 #mks
vexp = 10 * 1000 * 60 * 60 * 24 * 365 # m/yr






for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/CE34_69'):
    filenameC = files

cC = []
cvC = []
mC = []
mvC = []
dmC = []
dmvC = []
dnC = []
dnvC = []




for i in range(len(filenameC)):

    if True:

        wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/' + filenameC[i], usecols=0)
        flux = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/' + filenameC[i], usecols=1)
        v = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/' + filenameC[i], usecols=2)

	spl = splrep(wlen, flux, k = 3)
        f6 = splev(c6, spl)
        f9 = splev(c9, spl)

        splv = splrep(wlen, v, k = 3)
	v6 = splev(c6, splv)
        v9 = splev(c9, splv)
   

	int6 = simps(f6, c6)
	int9 = simps(f9, c9)

        int6up = simps(f6+v6, c6)
        int6low = simps(f6-v6, c6)
        int9up = simps(f9+v9, c9)
        int9low = simps(f9-v9, c9)


	int6 = int6/0.3
        int9 = int9/0.4

	int6up = int6up/0.3
        int9up = int9up/0.4

	int6low = int6low/0.3
        int9low = int9low/0.4
        
	
	mag6 = math.log(s6/int6, base)
        mag9 = math.log(s9/int9, base)
        color = mag6-mag9


	mag6up = math.log(s6/int6up, base)
        mag9up = math.log(s9/int9up, base)

	mag6low = math.log(s6/int6low, base)
        mag9low = math.log(s9/int9low, base)

	color_up = mag6up - mag9low
	color_low = mag6low - mag9up

        error = max(abs(color_up-color), abs(color_low-color))
	

        cC.append(color)
	cvC.append(error)
	
	
        gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
        
        tmpmlr = []
        tmpdustn = []
        tmpdustm = []
        for j in range(len(gram)):
	    mlr = g_data['MLR'][gram[j]] #Msun/yr
            rin = g_data['Rin'][gram[j]] #Rstar
            rout = rin * 1000 #Rstar
            mass = g_data['Mass'][gram[j]] * ms #kg
        
            gravity = (10**g_data['logg'][gram[j]])/100 #m/s2
     
            r = (G * mass / gravity)**0.5 * rout #m
            rsqr = G * mass / gravity * rin * rin #m2

            dustm = mlr  * r/vexp #g cm-2
            dustn = mlr *ms*1000 * r/vexp/4/np.pi/rsqr*1e-4/6e-15 #g cm-2


            tmpmlr.append(mlr)
            tmpdustm.append(dustm)
            tmpdustn.append(dustn)
        mvC.append(median_absolute_deviation(tmpmlr) * np.pi/(len(tmpmlr))**0.5)
        mC.append(tmpmlr[0])
	dnvC.append(median_absolute_deviation(tmpdustn) * np.pi/(len(tmpdustn))**0.5)
        dnC.append(tmpdustn[0])
	dmvC.append(median_absolute_deviation(tmpdustm) * np.pi/(30)**0.5)
        dmC.append(tmpdustm[0])


for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/ECAGB_69'):
    filenameE = files

cE = []
cvE = []
mE = []
mvE = []
dmE = []
dmvE = []
dnE = []
dnvE = []


i = 0
for i in range(len(filenameE)):

    if True:

        wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/ECAGB/' + filenameE[i], usecols=0)
        flux = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/ECAGB/' + filenameE[i], usecols=1)
        v = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/ECAGB/' + filenameE[i], usecols=2)

	spl = splrep(wlen, flux, k = 3)
        f6 = splev(c6, spl)
        f9 = splev(c9, spl)

        splv = splrep(wlen, v, k = 3)
	v6 = splev(c6, splv)
        v9 = splev(c9, splv)
   

	int6 = simps(f6, c6)
	int9 = simps(f9, c9)

        int6up = simps(f6+v6, c6)
        int6low = simps(f6-v6, c6)
        int9up = simps(f9+v9, c9)
        int9low = simps(f9-v9, c9)


	int6 = int6/0.3
        int9 = int9/0.4

	int6up = int6up/0.3
        int9up = int9up/0.4

	int6low = int6low/0.3
        int9low = int9low/0.4
        
	mag6 = math.log(s6/int6, base)
        mag9 = math.log(s9/int9, base)
        color = mag6-mag9


	mag6up = math.log(s6/int6up, base)
        mag9up = math.log(s9/int9up, base)

	mag6low = math.log(s6/int6low, base)
        mag9low = math.log(s9/int9low, base)

	color_up = mag6up - mag9low
	color_low = mag6low - mag9up

        error = max(abs(color_up-color), abs(color_low-color))
	

        cE.append(color)
	cvE.append(error)
        gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted(copy)/' + filenameE[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
        tmpmlr = []
        tmpdustn = []
        tmpdustm = []
        for j in range(len(gram)):
	    mlr = g_data['MLR'][gram[j]] #Msun/yr
            rin = g_data['Rin'][gram[j]] #Rstar
            rout = rin * 1000 #Rstar
            mass = g_data['Mass'][gram[j]] * ms #kg
        
            gravity = (10**g_data['logg'][gram[j]])/100 #m/s2
     
            r = (G * mass / gravity)**0.5 * rout #m
            rsqr = G * mass / gravity * rin * rin #m2

            dustm = mlr  * r/vexp #g cm-2
            dustn = mlr *ms*1000 * r/vexp/4/np.pi/rsqr*1e-4/6e-15 #g cm-2


            tmpmlr.append(mlr)
            tmpdustm.append(dustm)
            tmpdustn.append(dustn)
        mvE.append(median_absolute_deviation(tmpmlr) * np.pi/(len(tmpmlr))**0.5)
        mE.append(tmpmlr[0])
	dnvE.append(median_absolute_deviation(tmpdustn) * np.pi/(len(tmpdustn))**0.5)
        dnE.append(tmpdustn[0])
	dmvE.append(median_absolute_deviation(tmpdustm) * np.pi/(30)**0.5)
        dmE.append(tmpdustm[0])


for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/VRO_69'):
    filenameV = files

cV = []
cvV = []
mV = []
mvV = []
dmV = []
dmvV = []
dnV = []
dnvV = []


i = 0
for i in range(len(filenameV)):

    if True:

        wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/VRO/' + filenameV[i], usecols=0)
        flux = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/VRO/' + filenameV[i], usecols=1)
        v = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/VRO/' + filenameV[i], usecols=2)

	spl = splrep(wlen, flux, k = 3)
        f6 = splev(c6, spl)
        f9 = splev(c9, spl)

        splv = splrep(wlen, v, k = 3)
	v6 = splev(c6, splv)
        v9 = splev(c9, splv)
   

	int6 = simps(f6, c6)
	int9 = simps(f9, c9)

        int6up = simps(f6+v6, c6)
        int6low = simps(f6-v6, c6)
        int9up = simps(f9+v9, c9)
        int9low = simps(f9-v9, c9)


	int6 = int6/0.3
        int9 = int9/0.4

	int6up = int6up/0.3
        int9up = int9up/0.4

	int6low = int6low/0.3
        int9low = int9low/0.4
        
	mag6 = math.log(s6/int6, base)
        mag9 = math.log(s9/int9, base)
        color = mag6-mag9


	mag6up = math.log(s6/int6up, base)
        mag9up = math.log(s9/int9up, base)

	mag6low = math.log(s6/int6low, base)
        mag9low = math.log(s9/int9low, base)

	color_up = mag6up - mag9low
	color_low = mag6low - mag9up

        error = max(abs(color_up-color), abs(color_low-color))
	

        cV.append(color)
	cvV.append(error)
        gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
        tmpmlr = []
        tmpdustn = []
        tmpdustm = []
        for j in range(len(gram)):
	    mlr = g_data['MLR'][gram[j]] #Msun/yr
            rin = g_data['Rin'][gram[j]] #Rstar
            rout = rin * 1000 #Rstar
            mass = g_data['Mass'][gram[j]] * ms #kg
        
            gravity = (10**g_data['logg'][gram[j]])/100 #m/s2
     
            r = (G * mass / gravity)**0.5 * rout #m
            rsqr = G * mass / gravity * rin * rin #m2

            dustm = mlr  * r/vexp #g cm-2
            dustn = mlr *ms*1000 * r/vexp/4/np.pi/rsqr*1e-4/6e-15 #g cm-2


            tmpmlr.append(mlr)
            tmpdustm.append(dustm)
            tmpdustn.append(dustn)
        mvV.append(median_absolute_deviation(tmpmlr) * np.pi/(len(tmpmlr))**0.5)
        mV.append(tmpmlr[0])
	dnvV.append(median_absolute_deviation(tmpdustn) * np.pi/(len(tmpdustn))**0.5)
        dnV.append(tmpdustn[0])
	dmvV.append(median_absolute_deviation(tmpdustm) * np.pi/(len(tmpdustm))**0.5)
        dmV.append(tmpdustm[0])





axarr[0].errorbar(cE, mE, mvE, cvE, 'bo', ms = 20, lw = 3)
axarr[0].errorbar(cC, mC, mvC, cvC, 'ro', ms = 20, lw = 3)
axarr[0].errorbar(cV, mV, mvV, cvV, 'go', ms = 20, lw = 3)
axarr[0].legend(['CE0-2', 'CE3-4', 'VRO'], loc = 'best', fontsize = 35)
axarr[0].grid(True)
axarr[0].set_ylabel(r'dust mass loss rate [$M_{\odot}$yr$^{-1}$]')
axarr[0].set_yscale('log')

axarr[1].errorbar(cE, dnE, dnvE, cvE, 'bo', ms = 20, lw = 3)
axarr[1].errorbar(cC, dnC, dnvC, cvC, 'ro', ms = 20, lw = 3)
axarr[1].errorbar(cV, dnV, dnvV, cvV, 'go', ms = 20, lw = 3)
axarr[1].legend(['CE0-2', 'CE3-4', 'VRO'], loc = 'best', fontsize = 35)
axarr[1].grid(True)
axarr[1].set_ylabel(r'dust column density [cm$^{-2}$]')
axarr[1].set_yscale('log')

axarr[2].errorbar(cE, dmE, dmvE, cvE, 'bo', ms = 20, lw = 3)
axarr[2].errorbar(cC, dmC, dmvC, cvC, 'ro', ms = 20, lw = 3)
axarr[2].errorbar(cV, dmV, dmvV, cvV, 'go', ms = 20, lw = 3)
axarr[2].legend(['CE0-2', 'CE3-4', 'VRO'], loc = 'best', fontsize = 35)
axarr[2].grid(True)
axarr[2].set_ylabel(r'dust mass [$M_{\odot}$]')
axarr[2].set_yscale('log')


plt.savefig('3333.png')
plt.close()



        

