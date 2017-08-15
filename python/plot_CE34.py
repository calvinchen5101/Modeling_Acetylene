""" plot fits for all sources
to plot ECAGB and ERO, simply replace all CE34 in the script (cannot plot VRO), see VRO.py
"""


import ss_readfits
import numpy as np
from scipy.interpolate import splrep, splev
import scipy.constants as const
from astropy.stats import median_absolute_deviation
import matplotlib.pyplot as plt
import os

g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits")
tau_wlen = np.genfromtxt('/arrays/igloo1/ssp201701/relative_tau.dat', usecols = 0)
rel_tau = np.genfromtxt('/arrays/igloo1/ssp201701/relative_tau.dat', usecols = 1)

ms  = 1.98892e30 #kg
G = 6.67e-11
pc = const.parsec

loca = 10

def B(v, t):
    front = (2 * const.h * (v**3))/(const.c**2)
    den = np.e**(const.h * v/const.k/t) - 1
    #print np.e**(const.h * v/const.k/t)
    back = 1/den
    return front * back 

def IV0(v):
    return B(v, 3000)
    
def recal(i, iabs, v, t):
    #print ('i - bv')
    front_num = i - B(v, t)
    #print front_num
    back_num = iabs * IV0(v) - B(v, t)
    #print back_num
    den = IV0(v) - B(v, t)
    return (front_num * back_num / den + B(v, t)) * 1e26


for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted'):
    filename = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filename.append(files[i])


for i in range(len(filename)):
    wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/' + filename[i], usecols = 0)
    flux = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/' + filename[i], usecols = 1)
    v = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/' + filename[i], usecols = 2)
    allflux = []
    gram100 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filename[i], usecols = 0, skip_header = 1, max_rows = 30)
    c2h2100 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filename[i], usecols = 1, skip_header = 1, max_rows = 30)

    for j in range(len(gram100)):
	gwlen = g_data['Lspec'][gram100[j]]
	gflux = g_data['Fspec'][gram100[j]]
        cwlen = c_data['Lspec'][c2h2100[j]]
        cflux = c_data['Fspec'][c2h2100[j]]
        spl = splrep(gwlen, gflux, k=3)
        ggflux = splev(wlen, spl)
        spl = splrep(cwlen, cflux, k=3)
        ccflux = splev(wlen, spl)
        ctemp = c_data['T'][c2h2100[j]]

        gravity = (10**g_data['logg'][gram100[j]]) / 100 #m/s2
    
        mass = g_data['Mass'][gram100[j]] * ms
    
        rin = g_data['Rin'][gram100[j]]
        r_sqr = G*mass/gravity
        
        factor_out = (50000 * pc)**2/r_sqr/np.pi/loca/loca/rin/rin
        factor_photo = (50000 * pc)**2/r_sqr/np.pi
            
        photoflux = g_data['Fstar'][gram100[j]]
        tau = g_data['tau11_3'][gram100[j]]

        rel_tau = rel_tau * tau
        spl = splrep(tau_wlen, rel_tau, k = 3)
        rel_tau, obs_tau = splev(gwlen, spl), splev(wlen, spl)
    

        photoflux = photoflux * np.e**(-rel_tau *(loca-1)/loca)###
        gggflux = gflux * (1-np.e**(-rel_tau/loca))
        spl = splrep(gwlen, gggflux, k = 3)
        gggflux = splev(wlen, spl)
        gflux = gflux - photoflux
        spl = splrep(gwlen, gflux, k = 3)
        ngflux = splev(wlen, spl)
        spl = splrep(gwlen, photoflux, k = 3)
        npflux = splev(wlen, spl)
        intens_photo = npflux * 1e-26 * factor_photo
        intens_out = ngflux * 1e-26 * factor_out

        intens_sum = intens_photo + intens_out
        ratio = intens_photo / intens_sum
        wlen_expand = wlen * 1e-6
        
        niu = const.c/wlen_expand
        
        ncflux = recal(intens_sum, ccflux, niu, ctemp)
        comb_photo = ncflux * ratio / factor_photo 
        comb_out = ncflux * (1-ratio) / factor_out    
        comb = comb_photo + comb_out
        comb = comb * np.e**(-obs_tau/loca)
        comb = comb + gggflux
        allflux.append(comb)
    print filename[i]
    
    best = allflux[0]
    mad = median_absolute_deviation(allflux, axis = 0)
 
    fig = plt.figure(figsize = (10,7))
    ax1 = fig.add_subplot(1, 1, 1)
    ax2 = ax1.twinx()
    ax1.plot(wlen, flux, 'k-',lw = 3)
    ax1.plot(wlen, best, 'red',alpha = 0.75, lw = 3)
    ax2.plot(wlen, flux/best, 'orange', lw = 3)
    
    

    ax1.legend(['Observed', 'Best Fit'],loc = 'upper right')
    ax1.fill_between(wlen, flux - v, flux + v, facecolor = 'grey', alpha = 0.5, edgecolor = 'grey')
    ax1.fill_between(wlen, best - mad, best + mad, facecolor = 'pink', alpha = 0.5, edgecolor = 'pink')
    ax2.fill_between(wlen, flux/best - v/flux, flux/best + v/flux, facecolor = 'yellow', alpha = 0.25, edgecolor = 'yellow')
    ax2.axis([None,None,0,6])
    ax2.plot([0,40],[1,1],'k--', lw = 3)

    
   
    ax1.grid(True)

    plt.savefig('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filename[i][:-3] + 'png')
    #plt.show()
    plt.close()





        

    