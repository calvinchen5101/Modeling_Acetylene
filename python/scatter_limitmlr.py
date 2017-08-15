"""plot the limited MLR vs log N plot in the final presentation"""

import ss_readfits
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from astropy.stats import median_absolute_deviation
import matplotlib


g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits")
maxgas = np.genfromtxt('/arrays/igloo1/ssp201701/total_dust_m.txt', usecols = 2, skip_header = 1)

plt.xlabel('Column Density (log)')
plt.ylabel('Mass Loss Rate (log)')

matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize = (20,20))
plt.grid(True)
up = np.log10(1e1)
low = np.log10(1e1)

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/CE34_2012good'):
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

nallC = []
mallC = []
nall = []
mall = []

for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 0, skip_header = 1, max_rows = 30)
    
        
    tmpn = []
    tmpt = []
    tmpm = []
    jj = []
    for j in range(len(c2h2)):
	if c_data['logN'][c2h2[j]] <= maxgas[gram[j]] + up and c_data['logN'][c2h2[j]] >= maxgas[gram[j]] - low:
            tmpn.append(c_data['logN'][c2h2[j]])
            tmpt.append(c_data['T'][c2h2[j]])
            tmpm.append(np.log10(g_data['MLR'][gram[j]]))
	    jj.append(j)
    if len(tmpn)>0:
        nC.append(c_data['logN'][c2h2[jj[0]]])
        tC.append(c_data['T'][c2h2[jj[0]]])
        mC.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        nvC.append(median_absolute_deviation(tmpn))
        tvC.append(median_absolute_deviation(tmpt))
        mvC.append(median_absolute_deviation(tmpm))
	nallC.append(maxgas[gram[jj[0]]])
        mallC.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012good'):
    filenameE = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameE.append(files[i])
nE = []
tE = []
mE = []
nvE = []
tvE = [] 
mvE = [] 

nallE = []
mallE = []  

for i in range(len(filenameE)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 1, skip_header = 1, max_rows = 30)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 0, skip_header = 1, max_rows = 30)
    
    tmpn = []
    tmpt = []
    tmpm = []
    jj = []
    for j in range(len(c2h2)):
	if c_data['logN'][c2h2[j]] <= maxgas[gram[j]] + up and c_data['logN'][c2h2[j]] >= maxgas[gram[j]] - low:
            tmpn.append(c_data['logN'][c2h2[j]])
            tmpt.append(c_data['T'][c2h2[j]])
            tmpm.append(np.log10(g_data['MLR'][gram[j]]))
	    jj.append(j)
    if len(tmpn)>0:
        nE.append(c_data['logN'][c2h2[jj[0]]])
        tE.append(c_data['T'][c2h2[jj[0]]])
        mE.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        nvE.append(median_absolute_deviation(tmpn))
        tvE.append(median_absolute_deviation(tmpt))
        mvE.append(median_absolute_deviation(tmpm))
	nallE.append(maxgas[gram[jj[0]]])
        mallE.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/VRO_2012good'):
    filenameV = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameV.append(files[i])
nV = []
tV = []
mV = []
nvV = []
tvV = [] 
mvV = [] 

nallV = []
mallV = []  

for i in range(len(filenameV)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 1, skip_header = 1, max_rows = 30)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 0, skip_header = 1, max_rows = 30)
   
    tmpn = []
    tmpt = []
    tmpm = []
    jj = []
    for j in range(len(c2h2)):
	if c_data['logN'][c2h2[j]] <= maxgas[gram[j]] + up and c_data['logN'][c2h2[j]] >= maxgas[gram[j]] - low:
            tmpn.append(c_data['logN'][c2h2[j]])
            tmpt.append(c_data['T'][c2h2[j]])
            tmpm.append(np.log10(g_data['MLR'][gram[j]]))
	    jj.append(j)
    if len(tmpn)>0:
        nV.append(c_data['logN'][c2h2[jj[0]]])
        tV.append(c_data['T'][c2h2[jj[0]]])
        mV.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        nvV.append(median_absolute_deviation(tmpn))
        tvV.append(median_absolute_deviation(tmpt))
        mvV.append(median_absolute_deviation(tmpm))
	nallV.append(maxgas[gram[jj[0]]])
        mallV.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))



plt.plot(nall, mall, marker = '^', linestyle = 'none', color = 'grey', ms = 15)
plt.plot(nC, mC, marker = 'o', linestyle = 'none', color = 'red', ms = 15)
plt.plot(nE, mE, marker = 'o', linestyle = 'none', color = 'blue', ms = 15)

plt.plot(nV, mV, marker = 'o', linestyle = 'none', color = 'green', ms = 15)

for i in range(len(nC)):
     plt.plot([nC[i], nallC[i]], [mC[i], mallC[i]], linestyle = '-', marker = '.', color = '#ff9999', lw = 3)

for i in range(len(nE)):
     plt.plot([nE[i], nallE[i]], [mE[i], mallE[i]], linestyle = '-', marker = '.', color = '#99e6ff', lw = 3)

for i in range(len(nV)):
     plt.plot([nV[i], nallV[i]], [mV[i], mallV[i]], linestyle = '-', marker = '.', color = '#99ff99', lw = 3)


plt.xlabel('Column Density (log)')
plt.ylabel('Mass Loss Rate (log)')
plt.legend([ 'Derived', 'CE34', 'ECAGB', 'VRO'], loc = 'upper left')
plt.axis([17, 20.5, -9.5, -7])
plt.show()
    
    
