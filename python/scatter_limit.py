"""plot the gas T vs log N plot with a vertical line indicating the derived MLR in the final presentation
"""


import ss_readfits
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from astropy.stats import median_absolute_deviation


g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits")

maxgas = np.genfromtxt('/arrays/igloo1/ssp201701/total_dust_m.txt', usecols = 2, skip_header = 1)



for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/VRO_2012good'):
    filenameC = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameC.append(files[i])
nC = []
tC = []
mC = []
nvC = []
tvC = []    
mvC = []

for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameC[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
   
    tmpn = []
    tmpt = []

    tmpmaxgas = []
    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])
        tmpt.append(c_data['T'][c2h2[j]])
        
        tmpmaxgas.append(maxgas[gram[j]])
    plt.style.use('dark_background')
    plt.figure(figsize = (10,5))
    plt.plot([tmpmaxgas[0], tmpmaxgas[0]], [-100, 1000], color = 'orange', lw = 4)
    plt.plot(tmpn, tmpt, 'wo', ms = 8)
    plt.plot(tmpn[0], tmpt[0], 'magenta', marker = 'o', ms = 8 )
    
    plt.fill_between([min(tmpmaxgas), max(tmpmaxgas)], [-100, -100], [1000, 1000], facecolor = 'yellow', alpha = 0.25, edgecolor = 'yellow')
    
    plt.grid(True)
    plt.axis([14, 21, 0, 900])
    print filenameC[i]
    plt.show()
    

    
    
