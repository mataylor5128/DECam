'''Script for the calibration, reduction, and stacking
of multiple Dark Energy Camera images taken as part of the 
Survey of Centaurus A's Baryonic Structures observing campaign'''

#import numpy as np
#import matplotlib.pyplot as plt
#from astropy.io import fits, ascii
#from astropy import units as u
#from astropy.time import Time
#from astropy.coordinates import SkyCoord, EarthLocation, AltAz
#from astropy.stats import sigma_clip
import os
# from matplotlib.colors import LogNorm
# from pylab import cm
# import subprocess
#import sklearn
#from sklearn.linear_model import LinearRegression
# import pandas as pd
# from scipy.optimize import curve_fit
from photometric_calibration import photo_calib
from astrometric_calibration import astro_calib
from image_stacking          import image_stack
from background_subtract     import back_subtract
from create_mask             import make_mask

#import time

# exit()

# subprocess.call(['/Users/mtaylor/local/bin/./sex'])
# os.system('sex')
# subprocess.Popen('sex',shell=True)
# exit()

do_setup        = True
do_photo_calib  = False
do_astro_calib  = False
do_make_mask    = True
do_background   = True
do_stack_images = True
do_psfex        = True

#---------------------------#
#Setup files and directories#
#---------------------------#
prefix           = "/Users/mtaylor/Projects/SCABS"
calib_im_dir     = "%s/environ_study/raw_frames" % prefix
calib_out_dir    = "%s/pipeline_python/calib_out" % prefix
sci_im_dir       = "%s/CP_proc_sci_frames" % prefix
back_im_dir      = "%s/pipeline_python/back_subtracted_images" % prefix
SE_command       = "/Users/mtaylor/local/bin/./sex"
SE_config_dir    = "%s/pipeline_python/SE_config" % prefix
SE_out_dir       = "%s/pipeline_python/SE_output" % prefix
SCAMP_config_dir = "%s/pipeline_python/SCAMP_config" % prefix
SCAMP_out_dir    = "%s/pipeline_python/SCAMP_output" % prefix
BACK_out   = "%s/pipeline_python/BACKGROUND_output" % prefix
SWARP_config_dir = "%s/pipeline_python/SWARP_config" % prefix
SWARP_out_dir    = "%s/pipeline_python/SWARP_output" % prefix
PSFEx_config_dir = "%s/pipeline_python/PSFEx_config" % prefix
PSFEx_out_dir    = "%s/pipeline_python/PSFEx_output" % prefix
stack_im_dir     = "%s/pipeline_python/image_stacks" % prefix

filters          = ["i"]#,"g","r","i","z"]
tiles            = ["1"]

output_dirs = [calib_out_dir,SE_out_dir,SCAMP_out_dir,SWARP_out_dir,PSFEx_out_dir,back_im_dir,stack_im_dir,BACK_out]

if do_setup == True:
	#Create all the directories
	for dir in output_dirs:
		if not os.path.exists(dir): os.makedirs(dir)	
	
#----------------------------------------------------------#
#Compute Zeropoint, Airmass term, and Colour term Estimates#
#----------------------------------------------------------#
if do_photo_calib == True:
	photo_calib(calib_im_dir,calib_out_dir,sci_im_dir,SE_command,SE_config_dir,SE_out_dir,filters)

#-------------------------------------------------------------------#
#Perform the astrometric calibration using SCAMP with SE input files#
#-------------------------------------------------------------------#
if do_astro_calib == True:
	astro_calib(filters,tiles,SE_command,sci_im_dir,SE_config_dir,SE_out_dir,SCAMP_config_dir,SCAMP_out_dir,calib_out_dir)

#-----------------------------------------------------------------#
#Create source masks used for the background subtraction algorithm#
#-----------------------------------------------------------------#
if do_make_mask == True:
# 	do_tiles   = ["1","2","3","4","5","6","7"]
	do_tiles   = ["3"]
	do_filters = ["i"]
	do_dithers = ["1"]
	make_mask(do_tiles,do_filters,do_dithers,calib_out_dir,sci_im_dir,SE_config_dir,SE_out_dir,SE_command)
	exit()

#---------------------------------------------------------------------------#
#Perform dynamic spatial and time modelled background estimation/subtraction#
#---------------------------------------------------------------------------#
if do_background == True:
	back_type = "median"
	do_tile   = "1"
	do_filter = "i"
	do_dither = "1"
	n_backs   = 5
	back_subtract(back_type,do_tile,do_filter,do_dither,n_backs,sci_im_dir,SE_out_dir,BACK_out)
# 	tiles,SE_command,sci_im_dir,SE_config_dir,SE_out_dir,SCAMP_config_dir,SCAMP_out_dir,calib_out_dir)

#---------------------------------------------------------#
#Perform image stacking using SWARP with SCAMP input files#
#---------------------------------------------------------#
if do_stack_images == True:
	image_stack(filters,tiles,SWARP_config_dir,SWARP_out_dir,sci_im_dir,back_im_dir,stack_im_dir)

exit()

#Saturation = 65000 ADU

#SE for SCAMP input

#SCAMP for astrometric solution

#Background subtraction

#SWARP for final image stack