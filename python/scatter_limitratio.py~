"""plot MLR vs C in dust/C in gas with limited log N
"""


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
maxdust = np.genfromtxt('/asiaa/home/ssp201701/Downloads/dust_column.txt', usecols = 2, skip_header = 1)



matplotlib.rcParams.update({'font.size': 22})
plt.figure(figsize = (20,20))
plt.grid(True)
up = np.log10(1e9)
low = np.log10(1e9)

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/CE34_2012good'):
    filenameC = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameC.append(files[i])
nC = []
dC = []
mC = []
nvC = []
tvC = []    
mvC = []

nallC = []
mallC = []
nall = []
mall = []
dall = []

for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
    
    
        
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
        dC.append(maxdust[gram[jj[0]]])
        mC.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        nvC.append(median_absolute_deviation(tmpn))
        tvC.append(median_absolute_deviation(tmpt))
        mvC.append(median_absolute_deviation(tmpm))
	nallC.append(maxgas[gram[jj[0]]])
        mallC.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        dall.append(maxdust[gram[jj[0]]])
        

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012good'):
    filenameE = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameE.append(files[i])
nE = []
dE = []
mE = []
nvE = []
tvE = [] 
mvE = [] 

nallE = []
mallE = []  

for i in range(len(filenameE)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
    
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
        dE.append(maxdust[gram[jj[0]]])
        mE.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        nvE.append(median_absolute_deviation(tmpn))
        tvE.append(median_absolute_deviation(tmpt))
        mvE.append(median_absolute_deviation(tmpm))
	nallE.append(maxgas[gram[jj[0]]])
        mallE.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        dall.append(maxdust[gram[jj[0]]])

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/VRO_2012good'):
    filenameV = []
    for i in range(len(files)):
        if files[i][-4:] == '.txt':
            filenameV.append(files[i])
nV = []
dV = []
mV = []
nvV = []
tvV = [] 
mvV = [] 

nallV = []
mallV = []  

for i in range(len(filenameV)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 1, skip_header = 1, max_rows = 30, dtype = int)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 0, skip_header = 1, max_rows = 30, dtype = int)
    
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
        dV.append(maxdust[gram[jj[0]]])
        mV.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        nvV.append(median_absolute_deviation(tmpn))
        tvV.append(median_absolute_deviation(tmpt))
        mvV.append(median_absolute_deviation(tmpm))
	nallV.append(maxgas[gram[jj[0]]])
        mallV.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))
        dall.append(maxdust[gram[jj[0]]])

nC = np.array(nC)
dC = np.array(dC)
nC = nC - np.log10(6.02e23) + np.log10(26)
ratioC = dC - nC


nE = np.array(nE)
dE = np.array(dE)
nE = nE - np.log10(6.02e23) + np.log10(26)
ratioE = dE - nE

nV = np.array(nV)
dV = np.array(dV)
nV = nV - np.log10(6.02e23) + np.log10(26)
ratioV = dV - nV

nall = np.array(nall)
dall = np.array(dall)
nall = nall - np.log10(6.02e23) + np.log10(26)
ratioall = dall - nall






plt.errorbar(mC, ratioC, nvC, mvC, marker = 'o', linestyle = 'none', color = 'red', lw = 3, ms = 15)
plt.errorbar(mE, ratioE, nvE, mvE, marker = 'o', linestyle = 'none', color = 'blue', lw = 3, ms = 15)

plt.errorbar(mV, ratioV, nvV, mvV, marker = 'o', linestyle = 'none', color = 'green', lw = 3, ms = 15)
plt.plot(mall, ratioall, marker = 'o', linestyle = 'none', color = 'grey', ms = 15)

plt.xlabel('Mass Loss Rate (log)')
plt.ylabel('C in dust / C in C2H2 (log)')
plt.axis([None,None,-1.5,None])
plt.legend(['Derived', 'CE34', 'ECAGB', 'VRO'], loc = 'best')

plt.show()
    
    
