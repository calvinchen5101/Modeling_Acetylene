import ss_readfits
import numpy as np
from scipy.interpolate import splrep, splev
import os

g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012/10400.fits")
interwlen = wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/4001_spec.txt', usecols=0)


for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB_spectrum/ERO'):
    filename = files





iii = 0
while iii < len(filename):
    fff = open('/arrays/igloo1/ssp201701/CAGB/ERO_2012/' + filename[iii], 'a')

    
    
    wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/ERO/' + filename[iii], usecols=0)
    flux = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/ERO/' + filename[iii], usecols=1)
    v = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/ERO/' + filename[iii], usecols=2)
    
    if wlen[-1] < 20:
        clip = len(wlen)
    
    for l in range(len(wlen)):
        if wlen[l] > 20:
            clip = l
            break
    
 
    wlen = wlen[:clip]
    flux = flux[:clip]
    v = v[:clip]

    # loop through all dust model 
    rrr = []

    for i in range(len(g_data['GRAMS'])):
        gramflux = g_data['Fspec'][i]
        spl = splrep(interwlen, gramflux, k = 3)
        ngflux = splev(wlen, spl)
        
        r_lst = ((ngflux - flux)/v) ** 2
        r = r_lst.mean()
        rrr.append((str(g_data['GRAMS'][i]), str(g_data['C2H2'][i]), r))

    all_data = sorted(rrr, key=lambda r: r[2])
    for i in range(2000):
        fff.write(all_data[i][0] + '\t')
        fff.write(all_data[i][1] + '\t')
        fff.write('%.5E' % all_data[i][2] + '\n')
        #break ###
        print i
    fff.close()      
    print ('END')
    #break ###
    iii += 1


print ('REAL END')

#import sort_ERO
#import ERO_2000
    
        
    

        

    
