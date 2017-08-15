"""plot MLR vs C in dust/C in gas"""


import ss_readfits
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from astropy.stats import median_absolute_deviation


g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits")

maxgas = np.genfromtxt('/asiaa/home/ssp201701/Downloads/dust_column.txt', usecols = 2, skip_header = 1)



for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/CE34_2012good'):
    filenameC = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameC.append(files[i])
nC = []
dC = []
mC = []
nvC = []
  
mvC = []

for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
    nC.append(c_data['logN'][c2h2[0]])
    dC.append(maxgas[gram[0]])
    mC.append(np.log10(g_data['MLR'][gram[0]]))
    tmpn = []
   
    tmpm = []

    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])
       
        tmpm.append(np.log10(g_data['MLR'][gram[j]]))

   
    nvC.append(median_absolute_deviation(tmpn))

    mvC.append(median_absolute_deviation(tmpm))

nC = np.array(nC)
dC = np.array(dC)
nC = nC - np.log10(6.02e23) + np.log10(26)
ratio = dC - nC

plt.errorbar(mC, ratio, nvC, mvC, 'ro', ms = 10, lw = 3)

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012good'):
    filenameE = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameE.append(files[i])
nE = []
dE = []
mE = []
nvE = []
  
mvE = []

for i in range(len(filenameE)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
    nE.append(c_data['logN'][c2h2[0]])
    dE.append(maxgas[gram[0]])
    mE.append(np.log10(g_data['MLR'][gram[0]]))
    tmpn = []

    tmpm = []

    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])
        
        tmpm.append(np.log10(g_data['MLR'][gram[j]]))

    
    nvE.append(median_absolute_deviation(tmpn))

    mvE.append(median_absolute_deviation(tmpm))

nE = np.array(nE)
dE = np.array(dE)
nE = nE - np.log10(6.02e23) + np.log10(26)
ratio = dE - nE
plt.errorbar(mE, ratio, nvE, mvE, 'bo', ms = 10, lw = 3)

        

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/VRO_2012good'):
    filenameV = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameV.append(files[i])
nV = []
dV = []
mV = []
nvV = []
  
mvV = []

for i in range(len(filenameV)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
    nV.append(c_data['logN'][c2h2[0]])
    dV.append(maxgas[gram[0]])
    mV.append(np.log10(g_data['MLR'][gram[0]]))
    tmpn = []

    tmpm = []

    for j in range(len(c2h2)):
        tmpn.append(c_data['logN'][c2h2[j]])

        tmpm.append(np.log10(g_data['MLR'][gram[j]]))
 
    nvV.append(median_absolute_deviation(tmpn))

    mvV.append(median_absolute_deviation(tmpm))

nV = np.array(nV)
dV = np.array(dV)
nV = nV - np.log10(6.02e23) + np.log10(26)
ratio = dV - nV
plt.errorbar(mV, ratio, nvV, mvV, 'go', ms = 10, lw = 3)
plt.grid(True)
plt.legend(['CE34', 'ECAGB', 'VRO'], loc = 'best')
plt.xlabel('Dust Mass Lost Rate (log)')
plt.ylabel('C in dust / C in C2H2 (log)')
plt.show()


    
    
