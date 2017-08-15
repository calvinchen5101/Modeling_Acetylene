"""plot the saturation plot of C2H2 in the paper
"""


import ss_readfits
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
from matplotlib import rc


matplotlib.rcParams.update({'font.size': 35})
rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
data, tags = ss_readfits.ss_readfits('/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits')

allind = np.where(data['logN']==16)
allind = allind[0]
print allind


color = ['red', 'blue', 'green', 'orange']
t = []
plt.figure(figsize = (20,20))
for i in range(len(allind)):
    if data['T'][allind[i]] % 200 == 0:
        plt.plot(data['Lspec'][allind[i]], data['Fspec'][allind[i]], color = color[i/4], lw = 5, alpha = 0.75)
        t.append('T = ' + str(int(data['T'][allind[i]])) + ' K')

print data['T'][91]
plt.xlabel('wavelength [$\mu$m]')
plt.ylabel('normalized flux')
plt.grid(True)
plt.axis([6, 17, 0.4, 1.1])
plt.legend(t, loc = 'best')
plt.text(17, 0.4, r'N = 10$^{20}$ [cm$^{-2}$]', verticalalignment='bottom', horizontalalignment='right',color='black', fontsize=50)
plt.show()

