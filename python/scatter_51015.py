"""plot the 5, 10, 15 rin scatter plot in the paper
can be cahnged for the 2012 C2H2 model
"""


import ss_readfits
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from astropy.stats import median_absolute_deviation
import matplotlib
from matplotlib import rc

matplotlib.rcParams.update({'font.size': 35})

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/allconv_c2h2.fits")

plt.figure(figsize = (20,20))



filenameC = ['4244_spec.txt']
nC = []
tC = []

nvC = []
tvC = []    


for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_sorted5rin/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    
    tmpn = []
    tmpt = []
    
    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])
        tmpt.append(c_data['T'][c2h2[j]])
        
    plt.plot(tmpn, tmpt, 'rs', ms = 20, alpha = 0.5)
   
    plt.ylabel('gas temperature [K]')
    plt.xlabel(r'log (column density [cm$^{-2}$]')
    plt.grid(True)
    plt.axis([16, 20, 100, 800])
    print filenameC[i]

for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_sorted10rin/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    
    tmpn = []
    tmpt = []
    
    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])
        tmpt.append(c_data['T'][c2h2[j]])
        
    plt.plot(tmpn, tmpt, 'bv', ms = 20, alpha = 0.5)
    
    print filenameC[i]

for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_sorted15rin/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    
    tmpn = []
    tmpt = []
    
    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])
        tmpt.append(c_data['T'][c2h2[j]])
        
    plt.plot(tmpn, tmpt, 'g^', ms = 20,lw = 5, alpha = 0.5)
   

    print filenameC[i]
        


plt.text(17, 600, 'SSID: 4334\nGroup: CE3-4', verticalalignment='top', horizontalalignment='left',color='black', fontsize=35)

plt.legend([r'5 R$_{in}$', r'10 R$_{in}$', r'15 R$_{in}$'], loc = 'best', fontsize = 35)

plt.show()


    
    
