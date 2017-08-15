"""plot the gas T vs log N plot for all source
"""

import ss_readfits
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from astropy.stats import median_absolute_deviation


g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits")





for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012good'):
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
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30)
    
    tmpn = []
    tmpt = []

    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])
        tmpt.append(c_data['T'][c2h2[j]])
        
    plt.plot(tmpn, tmpt, 'ko')
    plt.plot(tmpn[0], tmpt[0], 'ro')
    plt.errorbar(np.median(tmpn), np.median(tmpt), median_absolute_deviation(tmpt), median_absolute_deviation(tmpn), 'go')
    
    plt.grid(True)
    plt.axis([14, 21, 0, 900])
    print filenameC[i]
    plt.savefig('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012good/' + filenameC[i][:-3] + 'png')
    plt.close()
    
        


    
    
