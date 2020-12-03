

import os
import argparse
import numpy as np
import astropy
from astropy.io import fits
import matplotlib
from matplotlib import pyplot as plt
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy import wcs
import math
import sys
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)



#Call in file.
mycube=sys.argv[1]
#Give me a list of coordinates in the form of RA, Dec in Hour Degrees.
cat_files = sys.argv[2]



#Call in the fits file and set the header and data columns.  Determine the shape to find which access is frequency.
#ASKAP data typically is frequency, stokes, dec, ra.  MWA data is stokes, frequency, dec, ra.
datacube = fits.open(mycube)
data = datacube[0].data
header = datacube[0].header
print data.shape


#Create a list of frequencies in the data cube and find the units of the frequency axis.
print datacube[0].header['CUNIT4']
rp = datacube[0].header['CRPIX4']
rf = datacube[0].header['CRVAL4']
df = datacube[0].header['CDELT4']
nf = datacube[0].header['NAXIS4']
xvals = rf + df*(np.arange(nf)-rp)
#xvals are the frequency in Hz np.subtract(xvals,1.66555e+09)
xvals=xvals[0:2999]
#print xvals



#Convert the frequency to velocity and make LSR Correction.  Get LSR sky frequency from 
#https://www.narrabri.atnf.csiro.au/observing/obstools/velo.html
#OH maser at rest frequency is 18cm or 0.18m wavelength
vels=np.multiply(np.subtract(1.665399e+09,xvals),0.180012)
#Correct for m/s into km/s
vels=np.divide(vels,1000)
#Reduce the number velocities to make the plot easier to read.
#vels=vels[500:999]
#print vels


#Open the text file with the ra and dec
file=open(cat_files)


#Assign the ra and dec variable names for each item in the list.
for line in file:
    ra,dec = line.strip().split(",")
    ra=str(ra)
    dec=str(dec)
    c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    print ra
    print dec



#convert the ra and dec to pixels on the image
    w = wcs.WCS(header, naxis=2)
    xpix,ypix=c.to_pixel(w,origin=0,mode='wcs')
    xpix=int(xpix)
    ypix=int(ypix)

    print xpix
    print ypix



#Build a list of values for the peak intensity at the pixel coordinates for each frequency.
    signal=[]
    for x in range(0, 3000):
        value = np.nanmean(data[x,:,ypix:ypix+1,xpix:xpix+1])
        #print value
        signal.append(value)
    

#Determine what is the maximum value (assuming the signal is in emission).
    max_signal=np.nan_to_num(signal)
    max_value = np.amax(max_signal)
    print max_value
#Determine the spectral RMS
    RMS=np.nanstd(signal[0:500])
    RMS2=str(round(RMS,2))
    print RMS
    
#Determine SNR for significance.
    snr=np.divide(max_value,RMS)
#Determine which channel holds the highest intensity peak.
    sliced=np.argmax(max_signal, axis=0)
    rms_number = np.nanstd(data[sliced,1000:1200,1000:1200])
    
    print sliced



    import matplotlib.patches as mpatches
#   Set the minor axis counters
    majorLocator = MultipleLocator(20)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(5)
    if np.isnan(signal[10])!=True:
    
#Make a spectra        
        bigfig=plt.figure(figsize=(20,12))
        ax1=bigfig.add_subplot(111)
        ax1.step(vels[0:2998],signal[0:2998],color='blue')
    #ax1.step(xvals[0:2000],signal[0:2000],color='red')
        ax1.set_title("SB10656"+" RA="+ra+" Dec="+dec, fontsize=24)
        ax1.set_xlabel("Velocity (km/s)",fontsize=24)
    #ax1.set_xlim(-60,30)
        ax1.set_ylabel("Flux Density (Jy)",fontsize=24)
        ax1.tick_params(labelsize=22, labelbottom=True)
        ax1.ticklabel_format(useOffset=False)
        ax1.xaxis.set_major_locator(majorLocator)
        ax1.xaxis.set_major_formatter(majorFormatter)
        plum_patch = mpatches.Patch(color='grey', label='1$\sigma$')
        ax1.fill_between(vels[0:3000], -0.022, 0.022, facecolor='grey', alpha=0.5)
        ax1.legend(handles=[plum_patch],  loc=1,fontsize=24)

# for the minor ticks, use no labels; default NullFormatter
        ax1.xaxis.set_minor_locator(minorLocator)

#Save the figure.     
        bigfig.savefig("OH Radical_"+str(xpix)+"_"+str(ypix)+".png")



