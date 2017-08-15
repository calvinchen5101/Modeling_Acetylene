from astropy.table import Table

def ss_readfits(fitsfile):
        if fitsfile[-5:] != '.fits': fitsfile+='.fits'
        data=Table.read(fitsfile)
        tags=data.colnames
        return data, tags