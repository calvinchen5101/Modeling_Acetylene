"""plot the fits for VRO
chane "fff" variable to the desire filename, can only do this one by one
"""


import os
import ss_readfits
import numpy as np
from astropy.stats import median_absolute_deviation
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy.interpolate import splrep, splev

plt.style.use('dark_background')

g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits")

#i_data, i_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/inform_c2h2.fits")
#data, tags = ss_readfits.ss_readfits("/data/ssp201701/Jan_combine2000.fits")
#print data['Fspec'][975*417+372]
for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted'):
    filename = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filename.append(files[i])

fff = '4246lo1_spec.txt'
plt.grid(True)
ms  = 1.98892e30 #kg
G = 6.67e-11
pc = const.parsec
rel_tau = np.genfromtxt('/arrays/igloo1/ssp201701/relative_tau.dat', usecols = 1)
tau_wlen = np.genfromtxt('/arrays/igloo1/ssp201701/relative_tau.dat', usecols = 0)
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

i = 0
while True:
    j = 0
    ttt = []
    #ind = np.where((c_data['T'] % 100 == 0) & (c_data['logN'] == 20))
    #ind = ind[0]
    ind = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + fff, usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    iii = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + fff, usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
    wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/VRO/' + fff, usecols = 0)
    flux = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/VRO/' + fff, usecols = 1)
    v = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/VRO/' + fff, usecols = 2)
    
    plt.plot(wlen, flux, 'w-')
    a = []
    print ind[0]
    while j < 30:
        gwlen = g_data['Lspec'][iii[j]]
        gflux = g_data['Fspec'][iii[j]]
        
        cwlen = c_data['Lspec'][ind[j]]
        cflux = c_data['Fspec'][ind[j]]
        spl = splrep(gwlen, gflux, k=3)
        ggflux = splev(wlen, spl)
        spl = splrep(cwlen, cflux, k=3)
        ccflux = splev(wlen, spl)
        ctemp = c_data['T'][ind[j]]

        gravity = (10**g_data['logg'][iii[j]]) / 100 #m/s2

        mass = g_data['Mass'][iii[j]] * ms
    
        rin = g_data['Rin'][iii[j]]
        r_sqr = G*mass/gravity
        
        factor_out = (50000 * pc)**2/r_sqr/np.pi/loca/loca/rin/rin
        factor_photo = (50000 * pc)**2/r_sqr/np.pi
            
        photoflux = g_data['Fstar'][iii[j]]
        tau = g_data['tau11_3'][iii[j]]

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
        
        a.append(comb)    
        ttt.append(ctemp)

    
       
        j+=1
    print fff
    best = a[0]
    mad = median_absolute_deviation(a, axis = 0)
    fig = plt.figure(figsize = (10,7))
    ax1 = fig.add_subplot(1, 1, 1)
    ax2 = ax1.twinx()
    ax1.plot(wlen, flux, 'w-', lw = 3)
    ax1.plot(wlen, best, 'magenta',alpha = 0.75, lw = 3)
    ax2.plot(wlen, flux/best, 'orange', lw = 3)
    
    
  
    ax1.legend(['Observed', 'Best Fit'],loc = 'best')
    ax1.fill_between(wlen, flux - v, flux + v, facecolor = 'grey', alpha = 0.5, edgecolor = 'grey')
    ax1.fill_between(wlen, best - mad, best + mad, facecolor = 'pink', alpha = 0.5, edgecolor = 'pink')
    ax2.fill_between(wlen, flux/best - v/flux, flux/best + v/flux, facecolor = 'yellow', alpha = 0.25, edgecolor = 'yellow')
    ax2.axis([None,None,0,6])
    ax2.plot([0,40],[1,1],'w--', lw = 3)

    
    
    ax1.grid(True)
   
    plt.show()
    break

