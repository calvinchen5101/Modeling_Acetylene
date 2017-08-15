import ss_readfits
import numpy as np
from scipy.interpolate import splrep, splev
import scipy.constants as const
from astropy.io import fits

wlen = np.genfromtxt('/arrays/igloo1/ssp201701/CAGB_spectrum/CE34/4001_spec.txt', usecols = 0)


g_data, g_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/grams_c.fits")
c_data, c_tags = ss_readfits.ss_readfits("/arrays/igloo1/ssp201701/fits/2012conv_C2H2.fits")

tau_wlen = np.genfromtxt('/arrays/igloo1/ssp201701/relative_tau.dat', usecols = 0)
rel_tau = np.genfromtxt('/arrays/igloo1/ssp201701/relative_tau.dat', usecols = 1)

ms  = 1.98892e30 #kg
G = 6.67e-11
pc = const.parsec

loca = 10

def B(v, t):
    front = (2 * const.h * (v**3))/(const.c**2)
    den = np.e**(const.h * v/const.k/t) - 1
    #print np.e**(const.h * v/const.k/t)
    back = 1/den
    return front * back 

def IV0(v):
    return B(v, 3000)
    
def recal(i, iabs, v, t):
    #print ('i - bv')
    front_num = i - B(v, t)
    #print front_num
    back_num = iabs * IV0(v) - B(v, t)
    #print back_num
    den = IV0(v) - B(v, t)
    return (front_num * back_num / den + B(v, t)) * 1e26

lenc = len(c_data['Lspec'])
lenw = len(wlen)
#leng = len(g_data['Teff'])

wlenc = c_data['Lspec'][0]


cflux = []
ctemp = c_data['T']
for i in range(len(c_data['Lspec'])):
    new_flux = c_data['Fspec'][i]
    spl = splrep(wlenc, new_flux, k=3)
    fflux = splev(wlen, spl)
    cflux.append(fflux)

cflux = np.array(cflux)
    
ctemp = ctemp + np.zeros((lenw, lenc))
ctemp = ctemp.transpose()
    
# loop through all dust model
gram = np.repeat(np.arange(10400, 11700), lenc)
c2h2 = np.tile(np.arange(lenc), 1300)
fspec = np.zeros((1, lenw))

for i in range(10400, 11700):
    gravity = (10**g_data['logg'][i]) / 100 #m/s2
    
    mass = g_data['Mass'][i] * ms
    
    rin = g_data['Rin'][i]
    r_sqr = G*mass/gravity
        
    factor_out = (50000 * pc)**2/r_sqr/np.pi/loca/loca/rin/rin
    factor_photo = (50000 * pc)**2/r_sqr/np.pi
    
        
    gflux = g_data['Fspec'][i]
    photoflux = g_data['Fstar'][i]
    tau = g_data['tau11_3'][i]

    rel_tau = rel_tau * tau
    spl = splrep(tau_wlen, rel_tau, k = 3)
    rel_tau, obs_tau = splev(g_data['Lspec'][i], spl), splev(wlen, spl)
    #print tau

    photoflux = photoflux * np.e**(-rel_tau *(loca-1)/loca)###
    spl = splrep(g_data['Lspec'][i], gflux, k = 3)
    gramsflux = splev(wlen, spl)
    gggflux = gflux * (1-np.e**(-rel_tau/loca))
    spl = splrep(g_data['Lspec'][i], gggflux, k = 3)
    gggflux = splev(wlen, spl)
    gflux = gflux - photoflux
    spl = splrep(g_data['Lspec'][i], gflux, k = 3)
    ngflux = splev(wlen, spl)
    spl = splrep(g_data['Lspec'][i], photoflux, k = 3)
    npflux = splev(wlen, spl)
    intens_photo = npflux * 1e-26 * factor_photo
    intens_out = ngflux * 1e-26 * factor_out

    intens_sum = intens_photo + intens_out
    ratio = intens_photo / intens_sum
    #ratio = 1 ###
    wlen_expand = wlen * 1e-6
        
    niu = const.c/wlen_expand
        
    ncflux = recal(intens_sum, cflux, niu, ctemp)
    comb_photo = ncflux * ratio / factor_photo 
    comb_out = ncflux * (1-ratio) / factor_out    
    comb = comb_photo + comb_out
    comb = comb * np.e**(-obs_tau/loca)
    comb = comb + gggflux
    #comb = comb_photo + ngflux ###
    fspec = np.concatenate((fspec, comb))
    print i



fspec = fspec[1:]


col1 = fits.Column(name='GRAMS', format = 'I', array=gram)
col2 = fits.Column(name='C2H2', format = 'I', array=c2h2)
col4 = fits.Column(name='Fspec', format = '365E', array=fspec)


cols = fits.ColDefs([col1, col2, col4])

tbhdu = fits.BinTableHDU.from_columns(cols)

tbhdu.writeto('/arrays/igloo1/ssp201701/fits/2012/11700.fits')
print ('end')


        

    
