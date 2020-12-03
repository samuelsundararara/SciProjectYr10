import astropy
from astropy.io import fits
from astropy.modeling import models
import numpy as np
import sys
import os.path
import matplotlib as mpl
mpl.use('Agg')
 
mpl.use()
mpl.use('Agg')
import matplotlib.pyplot as pyplot
import montage_wrapper as montage
#tables and votables
from astropy.io.votable import parse_single_table

from astropy.io import fits
from astropy.wcs import WCS

hdu = fits.open('G335_mean.fits')[0]
data = datacube[0].data
header = datacube[0].header
wcs = WCS(hdu.header)

print(hdu.shape)
new_cube=data[:,:,:,:].swapaxes(0,1)
rm_cube = np.squeeze(new_cube)
rm_cube.shape

datacube.writeto('G335_fixed.fits', clobber=True)
