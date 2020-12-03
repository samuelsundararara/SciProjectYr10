#!/bin/bash

#Call in particular add-ons that you want to use
import astropy
from astropy.io import fits
from astropy.modeling import models
import numpy as np
import matplotlib
from matplotlib import pyplot

#label the file you want to run the program on.
file1='image.restored.i.SB7797.cube.fits'

#open the data cube to do science on and set the data and header into different variables.
subcube=fits.open(file1)
data=subcube[0].data
header=subcube[0].header

#Set the cube that you want to work on and shape it.
rm_cube=data[:,:,:]

#set all zeros to nan to avoid it squeing the statistics.
rm_cube[rm_cube==0]=np.nan

#Run the mathmatical operation.  This one is taking the Standard Deviation.  The nanstd means it ignores all nans.
rms_cube=np.nanstd(rm_cube, axis=0)

#Create a signal to noise map.  Not applicable in all cases.  Toggle off using '#' in the front of the line.
snr=rm_cube/rms_cube

#Save the files of the map with the mathmatical operation and the new cube.  Change "G335_*" with appropriate name.  Make sure the file extension is *.fits.
subcube[0].data=np.float32(snr)
subcube.writeto('G335_snr.fits', clobber=True)
subcube[0].data=np.float32(rms_cube)
subcube.writeto('G335_rms.fits', clobber=True)
