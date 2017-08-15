"""create a convolved C2H2 model grid (IRS resolution)
"""

from astropy.io import fits
import numpy as np
from astropy.convolution import convolve, Gaussian1DKernel
import ss_readfits




c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012_C2H2.fits")

j = 0
while j < len(c_data['LAM'][0]):
    if c_data['Lspec'][0][j+1] > 5 and c_data['Lspec'][0][j] < 5:
        clip1 = j
    elif c_data['Lspec'][0][j] > 40:
        clip2 = j
	break
    j += 1


a1 = np.array(c_data['Lspec'][0][clip1:clip2]) + np.zeros((len(c_data['Lspec']), len(c_data['Lspec'][0][clip1:clip2])))
a2 = []
gauss = Gaussian1DKernel(5.5, x_size = 23)

for i in range(len(c_data['Fspec'])):
    smooth = convolve(c_data['Fspec'][i][clip1:clip2], gauss, boundary = 'extend', normalize_kernel = True)
    a2.append(smooth)
    print (i)


a2 = np.array(a2)
    
col1 = fits.Column(name='Lspec', format = '2496E', array=a1)
col2 = fits.Column(name='Fspec', format = '2496E', array=a2)
col3 = fits.Column(name='T', format = 'E', array=c_data['T'])
col4 = fits.Column(name='logN', format = 'E', array=c_data['LOGN'])
cols = fits.ColDefs([col3, col4, col1, col2])
print len(c_data['LAM'][0][clip1:clip2])
tbhdu = fits.BinTableHDU.from_columns(cols)

tbhdu.writeto('/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits')
print ('end')    

