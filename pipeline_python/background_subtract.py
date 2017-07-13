'''Script to perform a dynamic background estimation and subtraction
from multiple consecutive DECam images.'''

from astropy.io import fits
import numpy as np
import glob
from multiprocessing import Process, Queue
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from scipy import interpolate, stats
import itertools
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os,sys,subprocess
import time
#from random import randint
import random
# import inpaint
from email.mime.text import MIMEText
import smtplib
from astropy.stats import sigma_clip

def weight_im(weight_image,image):
	for chip_num in range(im_nchip):
		image[0][chip_num] = image[0][chip_num]*weight_image[0][chip_num]/np.max(weight_image[0][chip_num])
	return image

def find_nearest(array,nvals,val):
	sorted = np.sort(np.abs(array))
	keep_vals = sorted[0:nvals]
	inds = []
	for value in keep_vals:
		inds.append(list(np.abs(array)).index(value))
	inds = np.array(inds)
	return inds

def back_subtract(back_type,do_tile,do_filt,do_dither,n_backs,im_dir,SE_out_dir,back_dir):

	start_time = time.time()

	ss_im_out  = back_dir
	sky_im_out = back_dir
	mask_dir   = SE_out_dir

	#Identify the image to have background subtracted from
	im_file    = im_dir+"/sci_"+do_filt+"_p"+do_tile+"_d"+do_dither+"_image.fits"
	im_data    = fits.open(im_file)
	wt_file    = im_dir+"/sci_"+do_filt+"_p"+do_tile+"_d"+do_dither+"_wtmap.fits"
	wt_data    = fits.open(wt_file)
	mask_file  = mask_dir+"/sci_"+do_filt+"_p"+do_tile+"_d"+do_dither+"_image.SEG.fits"

	print "Removing problem chip..."
	#Remove chip S7 (extension 31) from calibration due to bad southern half of the chip
	im_frame = fits.HDUList()
	wt_frame = fits.HDUList()
	for ii in range(len(im_data)):
		if ii != 30:
	# 		print ii
			im_frame.append(im_data[ii])
			wt_frame.append(wt_data[ii])

	temp_image_file = "%s/temp_im_file.fits" % back_dir
	temp_wt_file    = "%s/temp_wt_file.fits" % back_dir

	im_frame.writeto(temp_image_file,overwrite=True)
	wt_frame.writeto(temp_wt_file,overwrite=True)
	im_data.close()
	wt_data.close()
	
	im_file = temp_image_file
	im_data = fits.open(temp_image_file)

	#Identify the set of frames to use for the background estimate
	#Determine the background files closest in time to the image
	if do_dither == '1':
		temp_back_ims = np.concatenate((glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+do_dither+"_image.fits"),glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)+1)+"_image.fits")))
		temp_back_weights = np.concatenate((glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+do_dither+'_wtmap.fits'),glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)+1)+"_wtmap.fits")))
		temp_mask_ims = np.concatenate((glob.glob(mask_dir+"/sci_"+do_filt+"_p*_d"+do_dither+"_image.SEG.fits"),glob.glob(mask_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)+1)+"_image.SEG.fits")))
	if do_dither != '1' and int(do_dither) <= 5:
		temp_back_ims = np.concatenate((glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)-1)+"_image.fits"),glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+do_dither+"_image.fits"),glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)+1)+"_image.fits")))
		temp_back_weights = np.concatenate((glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)-1)+"_wtmap.fits"),glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+do_dither+"_wtmap.fits"),glob.glob(im_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)+1)+"_wtmap.fits")))
		temp_mask_ims = np.concatenate((glob.glob(mask_dir+"/sci_"+do_filt+"_p*_d"+str(int(do_dither)-1)+"_image.SEG.fits"),glob.glob(mask_dir+"/sci_"+do_filt+"_p*_d"+do_dither+"_image.SEG.fits"),glob.glob(mask_dir+"sci_"+do_filt+"_t*_d"+str(int(do_dither)+1)+"_image.SEG.fits")))
	if int(do_dither) > 5:
		temp_back_ims = glob.glob(im_dir+"/sci_"+do_filt+"_p*_d*_image.fits")
		temp_back_weights = glob.glob(im_dir+"/sci_"+do_filt+"_p*_d*_wtmap.fits")
		temp_mask_ims = glob.glob(mask_dir+"/sci_"+do_filt+"_p*_d*_image.SEG.fits")

	im_h = fits.open(im_file)[0].header
	im_mjd = im_h['MJD-OBS']
	print "\nImage observation date = %f" % im_mjd

	#Delete tiles with bright, extended objects in them
	#SCABS: tile1 (CenA), tile11 (wCen), tile24(wCen)
	if do_dither == '1':
		gv_dithers = [do_dither,str(int(do_dither)+1)]
	if do_dither != '1' and int(do_dither) <= 5:
		gv_dithers = [str(int(do_dither)-1),do_dither,str(int(do_dither)+1)]
	if int(do_dither) > 5:
		gv_dithers = [str(ii+1) for ii in range(15)]

	bv_tiles = ['1',do_tile,'11','24']
	for gv_d in gv_dithers:
		for bv in bv_tiles:
			if im_dir+"/sci_"+do_filt+"_p"+bv+"_d"+gv_d+"_image.fits" in temp_back_ims:
				temp_back_ims = np.delete(temp_back_ims, list(temp_back_ims).index(im_dir+"/sci_"+do_filt+"_p"+bv+"_d"+gv_d+"_image.fits"),0)
				if im_dir+"/sci_"+do_filt+"_p"+bv+"_d"+gv_d+"_wtmap.fits" in temp_back_weights:
					temp_back_weights = np.delete(temp_back_weights, list(temp_back_weights).index(im_dir+"/sci_"+do_filt+"_p"+bv+"_d"+gv_d+"_wtmap.fits"), 0)
				if mask_dir+"/sci_"+do_filt+"_p"+bv+"_d"+gv_d+"_image.SEG.fits" in temp_mask_ims:
					temp_mask_ims = np.delete(temp_mask_ims, list(temp_mask_ims).index(mask_dir+"/sci_"+do_filt+"_p"+bv+"_d"+gv_d+"_image.SEG.fits"), 0)

	for imfile in temp_back_ims:
		if imfile.replace(im_dir,mask_dir).replace("image.fits","image.SEG.fits") not in temp_mask_ims:
			temp_back_ims = np.delete(temp_back_ims, list(temp_back_ims).index(imfile), 0)		
	for wfile in temp_back_weights:
		if wfile.replace(im_dir,mask_dir).replace("wtmap.fits","image.SEG.fits") not in temp_mask_ims:
			temp_back_weights = np.delete(temp_back_weights, list(temp_back_weights).index(wfile), 0)		

	back_mjds = []
	for ii in range(len(temp_back_ims)):
		back_mjds.append(fits.open(temp_back_ims[ii])[0].header["MJD-OBS"])
	back_mjds = np.array(back_mjds)

	time_diffs = back_mjds - im_mjd

	back_ims = []
	back_weights = []
	mask_ims = []

	for idx in find_nearest(time_diffs,n_backs,0):
		back_ims.append(temp_back_ims[idx])
		back_weights.append(temp_back_weights[idx])
		mask_ims.append(temp_mask_ims[idx])

	back_ims = np.array(back_ims)
	back_weights = np.array(back_weights)
	mask_ims = np.array(mask_ims)

	print "\nUsing the following frames for the background subtraction: " 
	print back_ims
	print "\nUsing the following frames for background masks: " 
	print mask_ims

# 	Create image to receive the final image
	im_final = fits.open(im_file)
	
	#Create image to receive the background map
	back_frame = fits.open(im_file)

	#Load 1st chip of image, , 5 "background" images, and 5 "masks"
	chip_nums = np.arange(1,len(im_data))#[1,2,3,4,5,6,7,8,9]
	for chip_num in chip_nums:
		print "Loading chip%i of image, backgrounds, and masks..." % chip_num
		im_chip       = im_data[chip_num].data
		back_chip_out = back_frame[chip_num].data 

		mask_chips = [fits.open(mask_ims[temp_mask])[chip_num].data for temp_mask in range(len(mask_ims))]
		if chip_num <= 29:
			back_chips = [fits.open(back_ims[temp_back])[chip_num].data for temp_back in range(len(back_ims))]
		else:
			back_chips = [fits.open(back_ims[temp_back])[chip_num+1].data for temp_back in range(len(back_ims))]

		#Iterate through the dpix x dpix image segments on all the background images
		#to calculate the medians and resulting means
		dpix = 256
		xchk = 0
		ychk = 0
		while ychk + dpix <= 4096:
			xchk = 0
			while xchk + dpix <= 2048:
# 				print "Section: [%i:%i,%i:%i]" % (xchk,xchk+dpix,ychk,ychk+dpix)
				temp_backs = [back_chips[ii][ychk:ychk+dpix,xchk:xchk+dpix] for ii in range(len(back_chips))]
				back_meds = [np.nanmedian(temp_backs[ii]) for ii in range(len(temp_backs))]
 				#Sigma-clip the medians to eliminate the worst points that escaped the masks
				back_meds = sigma_clip(back_meds,sigma=2.0)
				back_chip_out[ychk:ychk+dpix,xchk:xchk+dpix] = np.nanmean(back_meds)
				xchk += dpix
			ychk += dpix

		im_final[chip_num].data = im_chip - back_chip_out
		back_frame[chip_num].data = back_chip_out

	im_data.close()
	im_final.writeto("temp_back_subtracted.fits",overwrite=True)
	back_frame.writeto("temp_background.fits",overwrite=True)
	
	plt.figure(figsize=(10,10))
	for ii in range(16):#len(chip_nums)):
		plt.subplot(4,4,ii+1)
		plt.imshow(back_frame[chip_nums[ii]].data,origin='lower',cmap='viridis')
		plt.colorbar()
		
	plt.show()
	plt.close('all')
	exit()
	
	print "Plotting..."	
	plt.figure(1,figsize=(12,12)).clf()
	plt.imshow(im_chip,vmin=-2,vmax=5e3,origin='lower',cmap='viridis')

	plt.figure(2,figsize=(12,12)).clf()
	plt.subplot(3,2,1)
	plt.imshow(back_chips[0],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,2)
	plt.imshow(back_chips[1],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,3)
	plt.imshow(back_chips[2],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,4)
	plt.imshow(back_chips[3],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,5)
	plt.imshow(back_chips[4],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.title("Background Images")
	
	plt.figure(3,figsize=(12,12)).clf()
	plt.subplot(3,2,1)
	plt.imshow(mask_chips[0],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,2)
	plt.imshow(mask_chips[1],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,3)
	plt.imshow(mask_chips[2],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,4)
	plt.imshow(mask_chips[3],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.subplot(3,2,5)
	plt.imshow(mask_chips[4],vmin=-2,vmax=5e3,origin='lower',cmap='viridis')
	plt.title("Mask Images")
	plt.show()
	
	plt.close('all')
	exit()	
