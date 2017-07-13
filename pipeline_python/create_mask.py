'''Create a source mask of a DECam image using Source Extractor
segmentation maps and a convolution kernel'''

import glob
from astropy.io import fits, ascii
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.convolution import Gaussian2DKernel, convolve
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

def conv(image):
	kernel = Gaussian2DKernel(stddev=1)
	image = convolve(image,kernel)
	return image

def make_mask(do_tiles,do_filts,do_dithers,calib_out_dir,sci_im_dir,SE_config_dir,SE_out_dir,SE_command):
	for filter in do_filts:
		for do_tile in do_tiles:
			for do_dither in do_dithers:
				SE_config_file = SE_config_dir+"/ctio_decam_mask.sex"
				SE_param_file  = SE_config_dir+"/ctio_decam_mask.param"
				SE_filt_name   = SE_config_dir+"/default.conv"#gauss_3.5_7x7.conv"

				im_frameset = [sci_im_dir+"/sci_"+filter+"_p"+do_tile+"_d"+do_dither+"_image.fits"]
				wt_frameset = [ii.replace("image","wtmap") for ii in im_frameset]

				for ii in range(len(im_frameset)):
					print "Processing frame: %s" % im_frameset[ii]
			
					im_name = im_frameset[ii]
					wt_name = wt_frameset[ii]

					SE_cat_file = im_name.replace(sci_im_dir,SE_out_dir).replace(".fits","_mask.ldac")
					SE_xml_file = im_name.replace(sci_im_dir,SE_out_dir).replace(".fits","_mask.xml")
					SE_check_file = im_name.replace(sci_im_dir,SE_out_dir).replace(".fits",".BACK.fits")+","+im_name.replace(sci_im_dir,SE_out_dir).replace(".fits",".SEG.fits")

					#Open the image and weightmap
					temp_im_frame = fits.open(im_name)
					temp_wt_frame = fits.open(wt_name)

					#Get some basic header information
					im_h = temp_im_frame[0].header
					im_filter = im_h['FILTER'][0]
					im_exptime = im_h['EXPTIME']
					#Has the photometric calibration been performed? If so, load it.
					calib_files = glob.glob("%s/*calibration*.txt" % calib_out_dir)
					if len(calib_files) >= 1:
						calib_data = ascii.read(calib_files[-1])
						im_zp = calib_data['zp'][list(calib_data['filter']).index(filter)]
						im_k  = calib_data['k'][list(calib_data['filter']).index(filter)]
					else:
						im_zp = im_h['MAGZERO']
						im_k = 0.
					im_coord = SkyCoord(im_h['RA'],im_h['DEC'], unit=(u.hourangle, u.deg))
					im_ra = im_coord.ra.deg
					im_dec = im_coord.dec.deg
					im_mjd = im_h['MJD-OBS']

					#Calculate airmass of pointing
					ctio = EarthLocation(lat=im_h['OBS-LAT']*u.deg,lon=-1*im_h['OBS-LONG']*u.deg,height=im_h['OBS-ELEV']*u.m)
					time = Time(im_h['DATE-OBS'].replace('T',' '))
					im_altaz = im_coord.transform_to(AltAz(obstime=time,location=ctio,))
					im_airmass = im_altaz.secz

# 					kernel = Gaussian2DKernel(stddev=1)

					#Has SE produced the input source catalogue for this frame yet?
					#If not, extract the photometry
					if not os.path.exists(SE_cat_file):
						print "Removing problem chip..."
						#Remove chip S7 (extension 31) from calibration due to bad southern half of the chip
						im_frame = fits.HDUList()
						wt_frame = fits.HDUList()
# 						for ii in range(len(temp_im_frame)):
# 							if ii != 0:
# 								temp_im_frame[ii].data = convolve(temp_im_frame[ii].data,kernel)
# 								print "Convolving chip #%i" % ii
# 							if ii != 30:
# 						# 		print ii
# 								im_frame.append(temp_im_frame[ii])
# 								wt_frame.append(temp_wt_frame[ii])

						im_frame.append(temp_im_frame[0])
						wt_frame.append(temp_wt_frame[0])
						pool = Pool(processes=12)
						print "Convolving..."
						results = pool.map(conv, [temp_im_frame[ii].data for ii in range(1,61)])
						for ii in range(1,61):
							temp_im_frame[ii].data = results[ii-1]
							if ii != 30:
								im_frame.append(temp_im_frame[ii])
								wt_frame.append(temp_wt_frame[ii])

						#Check to ensure that multiprocessing hasn't changed the order of the chips
						for ii in range(1,61):
							if ii < 30:
# 								print ii, temp_im_frame[ii].header['EXTNAME'], im_frame[ii].header['EXTNAME']
								if temp_im_frame[ii].header['EXTNAME'] != im_frame[ii].header['EXTNAME']:
									print ii, temp_im_frame[ii].header['EXTNAME'], im_frame[ii].header['EXTNAME']
									print "Chips became unmatched!!"
									exit()
							if ii > 30:
# 								print ii, temp_im_frame[ii].header['EXTNAME'], im_frame[ii-1].header['EXTNAME']
								if temp_im_frame[ii].header['EXTNAME'] != im_frame[ii-1].header['EXTNAME']:
									print ii, temp_im_frame[ii].header['EXTNAME'], im_frame[ii-1].header['EXTNAME']
									print "Chips became unmatched!!"
									exit()

						im_frame.writeto('temp_im_file.fits',overwrite=True)
						wt_frame.writeto('temp_wt_file.fits',overwrite=True)

						#Get the median FWHM from the chipset
						temp_fwhms = []
						for ii in range(len(im_frame)):
							if ii != 0:
								temp_fwhms.append(im_frame[ii].header['FWHM'])
						im_fwhm = np.median(temp_fwhms)*0.263

						#Last SE details, create and execute SE command
						SE_satur_level = 45000
						im_gain = 4.0
						im_rdnoise = 6.0

						print "Extracting sources..."
						command = [SE_command,'temp_im_file.fits','-c',SE_config_file,'-CATALOG_NAME',SE_cat_file,'-PARAMETERS_NAME',SE_param_file,'-WEIGHT_IMAGE','temp_wt_file.fits','-FILTER_NAME',SE_filt_name,'-XML_NAME',SE_xml_file,'-SATUR_LEVEL',str(SE_satur_level),'-MAG_ZEROPOINT',str(im_zp),'-SEEING_FWHM',str(im_fwhm),'-CHECKIMAGE_TYPE','BACKGROUND,SEGMENTATION','-CHECKIMAGE_NAME',SE_check_file]
						subprocess.call(command)