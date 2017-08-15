import os
import numpy as np

for root, dirs, files in os.walk('/arrays/igloo1/ssp201701/CAGB/CE34_2012'):
    filename = files
    
iii = 0
while iii < len(filename):
    gram = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012/' + filename[iii], usecols=0, skip_header = 1, dtype = str)
    c2h2 = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012/' + filename[iii], usecols=1, skip_header = 1, dtype = str)
    chi = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB/CE34_2012/' + filename[iii], usecols=2, skip_header = 1)
    
    r_tuples = []
    for i in range(len(gram)):
        r_tuples.append((gram[i], c2h2[i], chi[i]))

    all_data = sorted(r_tuples, key=lambda r: r[2])

    fff = open('/arrays/igloo1/ssp201701/CAGB/CE34_2012sorted/' + filename[iii], 'w')
    fff.write('GRAM\t')
    fff.write('c2h2\t')
    fff.write('chi_tot\n')

    for i in range(2000):
        fff.write(str(all_data[i][0]) + '\t')
        fff.write(str(all_data[i][1]) + '\t')
        fff.write('%.5E' % all_data[i][2] + '\n')
        #fff.write(all_data[i][3] + '\n')
    
    fff.close()
    print (iii)
    iii += 1

    
    
