#!/bin/bash

#Call in your required add-ons
import astropy
from astropy.io import fits
from astropy.modeling import models
import numpy as np
import sys
import os.path

#Call-in named files on command line
file1=sys.argv[1]

#open up the fits file and set the data and header columns
datacube=fits.open(file1)
data=datacube[0].data
header=datacube[0].header
#set the data shape and axis
cube=data[:,:,:,:]

#calculate the standard deviation of channel 10 at the pixel position of 400:700
rms_number=np.nanstd(cube[10:11,1:1,400:700,400:700])

#print the value to screen
print rms_number

#save the value to a file.  If the file doesn't exist, create the file first.
if os.path.exists('rms_final'):
    with open('rms_final','a') as f:
        f.write(str(rms_number)+"\n")
else:
    f=open('rms_final','w')
    f.write(str(rms_number)+"\n")


