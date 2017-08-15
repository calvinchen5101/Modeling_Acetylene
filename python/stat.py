"""generate txt file with parameters
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
plt.style.use('dark_background')


plt.figure(figsize = (18,20))
plt.grid(True)
up = np.log10(1e4)
low = np.log10(1e4)
fff = open('stat_exp.txt', 'w')
fff.write('SSID\tMLR\tlogN\tMLR_err\tlogN_err\n')
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

nall = []
mall = []

for i in range(len(filenameC)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 1, skip_header = 1, max_rows = 30)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filenameC[i], usecols = 0, skip_header = 1, max_rows = 30)
    #t1 = []
    #t2 = []
    #c2h2 = list(c2h2)
    #for j in range(len(c2h2)):
	#if c2h2.count(c2h2[j]) >= 3:
	    #t1.append(c2h2[j])
    	    #t2.append(gram[j])
    #c2h2 = t1
    #gram = t2
    
        
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
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	fff.write(filenameC[i][:-9] + '\t')
	fff.write(str(round(tmpm[0],3)) + '\t')
        fff.write(str(round(tmpn[0], 1)) + '\t')
        fff.write(str(round(median_absolute_deviation(tmpm),3)) + '\t')
        fff.write(str(round(median_absolute_deviation(tmpn),3)) + '\n')

        

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

for i in range(len(filenameE)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 1, skip_header = 1, max_rows = 30)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/ECAGB_2012sorted/' + filenameE[i], usecols = 0, skip_header = 1, max_rows = 30)
    #t1 = []
    #t2 = []
    #c2h2 = list(c2h2)
    #for j in range(len(c2h2)):
	#if c2h2.count(c2h2[j]) >= 3:
	    #t1.append(c2h2[j])
    	    #t2.append(gram[j])
    #c2h2 = t1
    #gram = t2
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
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	fff.write(filenameE[i][:-9] + '\t')
	fff.write(str(round(tmpm[0],3)) + '\t')
        fff.write(str(round(tmpn[0], 1)) + '\t')
        fff.write(str(round(median_absolute_deviation(tmpm),3)) + '\t')
        fff.write(str(round(median_absolute_deviation(tmpn),3)) + '\n')

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

for i in range(len(filenameV)):
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 1, skip_header = 1, max_rows = 30)
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/VRO_2012sorted/' + filenameV[i], usecols = 0, skip_header = 1, max_rows = 30)
    #t1 = []
    #t2 = []
    #c2h2 = list(c2h2)
    #for j in range(len(c2h2)):
	#if c2h2.count(c2h2[j]) >= 3:
	    #t1.append(c2h2[j])
    	    #t2.append(gram[j])
    #c2h2 = t1
    #gram = t2
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
	nall.append(maxgas[gram[jj[0]]])
        mall.append(np.log10(g_data['MLR'][gram[jj[0]]]))
	fff.write(filenameV[i][:-9] + '\t')
	fff.write(str(round(tmpm[0],3)) + '\t')
        fff.write(str(round(tmpn[0], 1)) + '\t')
        fff.write(str(round(median_absolute_deviation(tmpm),3)) + '\t')
        fff.write(str(round(median_absolute_deviation(tmpn),3)) + '\n')

#plt.errorbar(nall, mall, None, None, marker = 'o', linestyle = 'none', color = 'grey', lw = 2, ms = 10)
plt.errorbar(nC, tC, tvC, nvC, marker = 'o', linestyle = 'none', color = '#ff9999', lw = 2, ms = 10)
plt.errorbar(nE, tE, tvE, nvE, marker = 'o', linestyle = 'none', color = '#99e6ff', lw = 2, ms = 10)
#plt.errorbar(0,0,0,0,marker = 'o', linestyle = 'none', color = '#99e6ff', lw = 2, ms = 10)
plt.errorbar(nV, tV, tvV, nvV, marker = 'o', linestyle = 'none', color = '#99ff99', lw = 2, ms = 10)

plt.legend([ 'CE34', 'ECAGB', 'VRO'], loc = 'upper left')
plt.axis([17, 20.5, 0, 900])
#plt.show()
fff.close()
    
    
