'''Use SCAMP and Source Extractor to estimate
the astrometric calibration for a set of DECam
science images.'''

def astro_calib(filters,SE_command,sci_im_dir,SE_config_dir,SE_out_dir,SCAMP_config_dir,SCAMP_out_dir,calib_out_dir):
	import glob
	import numpy as np
	from astropy.io import fits, ascii
	from astropy import units as u
	from astropy.time import Time
	from astropy.coordinates import SkyCoord, EarthLocation, AltAz
	import os
	import subprocess
	import select_stars as ss
	import matplotlib.pyplot as plt
	import re

# 	print os.environ["DYLD_LIBRARY_PATH"]
# 	subprocess.call(["convert"])
# 	exit()
	
	SE_config_file = SE_config_dir+"/ctio_decam_scamp.sex"
	SE_param_file  = SE_config_dir+"/ctio_decam_scamp.param"
	SE_filt_name   = SE_config_dir+"/default.conv"#gauss_3.5_7x7.conv"

	#Set up the SCAMP parameters
	scamp_pos_err       = "1.2"
	scamp_scale_err     = "1.02"
	scamp_angle_err     = "0.02"
	scamp_sn_thresh     = "40.,80."
	scamp_fwhm_thresh   = "0.,100."
	scamp_crossid_rad   = "4."
	scamp_distort_deg   = "4"
	scamp_photclip_nsig = "2."
	scamp_astref_cat    = "2MASS"
	scamp_astref_band   = "Ks"
	scamp_astref_maglim = "6.,20."
	scamp_match_resol   = "0."
	
	scamp_refcat_dir    = "%s/refcat" % SCAMP_out_dir
	#Create the reference catalogue directory if necessary
	if not os.path.exists(scamp_refcat_dir):
		os.makedirs(scamp_refcat_dir)
	scamp_cat_type_out  = "FITS_LDAC"
	scamp_check_type    = "SKY_ALL,FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR"
	scamp_check_file    = ",".join([SCAMP_out_dir+"/"+mystring for mystring in ["sky_all","fgroups","distort","astr_interror2d","astr_interror1d","astr_referror2d","astr_referror1d","astr_chi2","psphot_error"]])

	#Loop through each requested filter
	for filter in filters:
		im_frameset = glob.glob("%s/sci_%s_*image.fits" % (sci_im_dir,filter))
# 		im_frameset = [im_frameset[0],im_frameset[1]]
		wt_frameset = [ii.replace("image","wtmap") for ii in im_frameset]
		
		scamp_output_dir    = "%s/PROG2014A-0610_FILTER%s_ORDER%s_REFCAT%s" % (SCAMP_out_dir,filter,scamp_distort_deg,scamp_astref_cat)
		scamp_list_file     = "%s/scamp_decam_2014A-0610_%s.dat" % (SCAMP_out_dir,filter)
		scamp_cat_file_out  = "%s/scamp_decam_2014A-0610_%s.ldac" % (SCAMP_out_dir,filter)
		scamp_list_file     = "%s/2014A-0610_FILTER%s_images.SCAMP.list" % (SCAMP_config_dir,filter)

		for ii in range(len(im_frameset)):
			print "Processing frame: %s" % im_frameset[ii]
			
			im_name = im_frameset[ii]
			wt_name = wt_frameset[ii]

			SE_cat_file = im_name.replace(sci_im_dir,SE_out_dir).replace(".fits","_scamp.ldac")
			SE_xml_file = im_name.replace(sci_im_dir,SE_out_dir).replace(".fits","_scamp.xml")
			SE_check_file = im_name.replace(sci_im_dir,SE_out_dir).replace(".fits",".BACK.fits")

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

			#Has SE produced the input source catalogue for this frame yet?
			#If not, extract the photometry
			if not os.path.exists(SE_cat_file):
				print "Removing problem chip..."
				#Remove chip S7 (extension 31) from calibration due to bad southern half of the chip
				im_frame = fits.HDUList()
				wt_frame = fits.HDUList()
				for ii in range(len(temp_im_frame)-1):
					if ii != 30:
				# 		print ii
						im_frame.append(temp_im_frame[ii])
						wt_frame.append(temp_wt_frame[ii])
				im_frame.writeto('temp_im_file.fits',clobber=True)
				wt_frame.writeto('temp_wt_file.fits',clobber=True)

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
				command = [SE_command,'temp_im_file.fits','-c',SE_config_file,'-CATALOG_NAME',SE_cat_file,'-PARAMETERS_NAME',SE_param_file,'-WEIGHT_IMAGE','temp_wt_file.fits','-FILTER_NAME',SE_filt_name,'-XML_NAME',SE_xml_file,'-SATUR_LEVEL',str(SE_satur_level),'-MAG_ZEROPOINT',str(im_zp),'-SEEING_FWHM',str(im_fwhm),'-CHECKIMAGE_TYPE','BACKGROUND','-CHECKIMAGE_NAME',SE_check_file]
				subprocess.call(command)

			#Free up memory
			temp_im_frame.close()
			temp_wt_frame.close()

			#Set up the SCAMP .head and .ahead files
			scamp_cat_file   = "%s/2014A-0610_FILTER%s" % (SCAMP_out_dir,im_name.replace(sci_im_dir+"/sci_","").replace("image.fits","stars.ldac"))
			scamp_head_file  = scamp_cat_file.replace(".ldac",".head")
			scamp_ahead_file = scamp_head_file.replace(".head",".ahead")

			#Does the SCAMP star catalogue already exist for this image?
			print "Looking for or creating star catalogue: %s" % scamp_cat_file
			if not os.path.exists(scamp_cat_file):
				#Open SE catalogue and create a new empty catalogue.
				temp_SE_cat = fits.open(SE_cat_file,ignore_missing_end=True)
				temp_star_cat = fits.HDUList()
				for jj in range(len(temp_SE_cat)):
					if jj == 0:
						temp_star_cat.append(temp_SE_cat[jj])
						temp_star_cat[jj].header = temp_SE_cat[jj].header
					else:
						if jj%2 != 0: #is jj odd?
							temp_star_cat.append(temp_SE_cat[jj])
							temp_star_cat[jj].header = temp_SE_cat[jj].header
# 							temp_star_cat[jj].header.set('CDELT1','0.0000731')
# 							temp_star_cat[jj].header.set('CUNIT1','deg')
# 							temp_star_cat[jj].header.set('CTYPE1','RA---TAN')
# 							temp_star_cat[jj].header.set('CDELT2','0.0000731')
# 							temp_star_cat[jj].header.set('CUNIT2','deg')
# 							temp_star_cat[jj].header.set('CTYPE2','DEC--TAN')
# 							temp_star_cat[jj].header.set('EQUINOX','2000.0')
# 							temp_star_cat[jj].header.set('INSTRUME','DECam')
						else:
							#Only choose sources with FLAG = 0
							temp_data = temp_SE_cat[jj].data[(temp_SE_cat[jj].data['FLAGS'] == 0)]
							gv_stars = ss.select_stars(temp_data['FLUX_RADIUS'],temp_data['MAG_AUTO'])
							gv_star_mask = np.zeros(len(temp_data),dtype=bool)
							for ind in gv_stars:
								gv_star_mask[ind] = True
							temp_data = temp_data[gv_star_mask]
							if len(temp_data) < 10:
								print "Not enough stars on chip %s!" % temp_SE_cat[jj].header['EXTNAME']						
							temp_hdu = fits.BinTableHDU.from_columns(temp_SE_cat[jj].columns,nrows=len(temp_data))
							temp_hdu.header = temp_SE_cat[jj].header
							temp_hdu.data = temp_data
							temp_star_cat.append(temp_hdu)
	# 							temp_star_cat[jj].header.set('CTYPE1','RA---TAN')
	# 							temp_star_cat[jj].header.set('CTYPE2','RA---TAN')
	# 							temp_star_cat[jj].header.set('EQUINOX','2000.0')
				temp_star_cat.writeto(scamp_cat_file,clobber=True)		

# 			print "Looking for or creating star catalogue: %s" % scamp_cat_file
# 			if not os.path.exists(scamp_cat_file):
# 				#Open SE catalogue and create a new empty catalogue.
# 				temp_SE_cat = fits.open(SE_cat_file)
# 				temp_star_cat = fits.HDUList()
# 				for jj in [x for x in np.arange(2,117) if x%2 ==0]:
# # 					if jj == 0:
# # 						temp_star_cat.append(temp_SE_cat[jj])
# # 					else:
# # 						if jj%2 != 0: #is jj odd?
# # 							temp_star_cat.append(temp_SE_cat[jj])
# # 						else:
# 					#Only choose sources with FLAG = 0
# 					temp_data = temp_SE_cat[jj].data[(temp_SE_cat[jj].data['FLAGS'] == 0)]
# 					gv_stars = ss.select_stars(temp_data['FLUX_RADIUS'],temp_data['MAG_AUTO'])
# 					gv_star_mask = np.zeros(len(temp_data),dtype=bool)
# 					for ind in gv_stars:
# 						gv_star_mask[ind] = True
# 					temp_data = temp_data[gv_star_mask]
# 					if len(temp_data) < 10:
# 						print "Not enough stars on chip %s!" % temp_SE_cat[jj].header['EXTNAME']						
# 					temp_hdu = fits.BinTableHDU.from_columns(temp_SE_cat[jj].columns,nrows=len(temp_data))
# 					temp_hdu.header = temp_SE_cat[jj].header
# 					temp_hdu.data = temp_data
# 					temp_star_cat.append(temp_hdu)
# 				temp_star_cat.writeto(scamp_cat_file,clobber=True)		

			#Does the SCAMP .ahead already exist for this image?
			print "Looking for or creating SCAMP .ahead file: %s" % scamp_ahead_file
			if not os.path.exists(scamp_ahead_file):
				temp_scamp_ahead_file = open(scamp_ahead_file,'w')
				for jj in [x for x in np.arange(2,117) if x%2 ==0]:
				
					print >> temp_scamp_ahead_file, "PROGRAM = '2014A-0610'"
					print >> temp_scamp_ahead_file, "TYPE    = 'science'"
					print >> temp_scamp_ahead_file, "FILTER  = '%s'" % filter
					print >> temp_scamp_ahead_file, "AIRMASS = %.2f" % im_airmass
					print >> temp_scamp_ahead_file, "EXPTIME = %i" % im_exptime
					print >> temp_scamp_ahead_file, "PHOTFLAG= 'F'"
					print >> temp_scamp_ahead_file, "PHOT_ZP = %.2f" % im_zp
					print >> temp_scamp_ahead_file, "PHOT_K  = %.2f" % im_k
					print >> temp_scamp_ahead_file, "CUNIT1  = 'deg'"
					print >> temp_scamp_ahead_file, "CTYPE1  = 'RA---TAN'" 
					print >> temp_scamp_ahead_file, "CTYPE2  = 'deg'"
					print >> temp_scamp_ahead_file, "CTYPE2  = 'DEC--TAN'"
					print >> temp_scamp_ahead_file, "EQUINOX = 2000.0"
					print >> temp_scamp_ahead_file, "INSTRUME= 'DECam'"
					print >> temp_scamp_ahead_file, "END"
				temp_scamp_ahead_file.close()
				
			#Does the image list file exist? If not, create it. Otherwise add to it.
			print "Looking for or creating SCAMP image list: %s" % scamp_list_file
			if not os.path.exists(scamp_list_file):
				temp_scamp_list_file = open(scamp_list_file,'a')
				print >> temp_scamp_list_file, scamp_cat_file
			else:
				if "%s\n" % scamp_cat_file not in file(scamp_list_file).readlines():
					temp_scamp_list_file = open(scamp_list_file,'a')
					print >> temp_scamp_list_file, scamp_cat_file
					temp_scamp_list_file.close()

# 		command='scamp @'+scamp_list_file+' -c scamp_config/ctio_decam.scamp'+' -MERGEDOUTCAT_TYPE '+scamp_cat_type_out+' -MERGEDOUTCAT_NAME '+scamp_cat_file_out+' -MATCH Y -WRITE_XML Y -XML_NAME '+scamp_xml_file+' -SAVE_REFCATALOG Y -REFOUT_CATPATH '+scamp_refcat_dir+' -CHECKPLOT_DEV PSC -CHECKPLOT_ANTIALIAS Y -CHECKPLOT_TYPE '+scamp_check_type + ' -CHECKPLOT_NAME '+scamp_check_file + ' -ASTREF_CATALOG '+scamp_astref_catalog+' -ASTREF_BAND '+scamp_astref_band+' -ASTREFMAG_LIMITS '+scamp_astrefmag_limits+' -DISTORT_DEGREES '+scamp_distort_degrees+' -PHOTCLIP_NSIGMA 2. -SOLVE_ASTROM Y -SOLVE_PHOTOM Y -POSITION_MAXERR '+scamp_pos_error+' -PIXSCALE_MAXERR '+scamp_scale_error+' -POSANGLE_MAXERR '+scamp_angle_error+' -SN_THRESHOLDS '+scamp_sn_thresholds+' -FWHM_THRESHOLDS '+scamp_fwhm_thresholds+' -CROSSID_RADIUS '+scamp_crossid_radius+' -MATCH_RESOL '+scamp_match_resol
		command = ['/Users/mtaylor/local/bin/scamp','@'+scamp_list_file,'-c',SCAMP_config_dir+"/ctio_decam.scamp","-MERGEDOUTCAT_TYPE",scamp_cat_type_out,"-MERGEDOUTCAT_NAME",scamp_cat_file_out,"-MATCH","Y","-WRITE_XML","Y","-XML_NAME","scamp.xml","-SAVE_REFCATALOG","Y","-REFOUT_CATPATH",scamp_refcat_dir,"-CHECKPLOT_DEV","PSC","-CHECKPLOT_ANTIALIAS","Y","-CHECKPLOT_TYPE",scamp_check_type,"-CHECKPLOT_NAME",scamp_check_file,"-ASTREF_CATALOG",scamp_astref_cat,"-ASTREF_BAND",scamp_astref_band,"-ASTREFMAG_LIMITS",scamp_astref_maglim,"-DISTORT_DEGREES",scamp_distort_deg,"-PHOTCLIP_NSIGMA",scamp_photclip_nsig,"-SOLVE_ASTROM","Y","-POSITION_MAXERR",scamp_pos_err,"-PIXSCALE_MAXERR",scamp_scale_err,"-POSANGLE_MAXERR",scamp_angle_err,"-SN_THRESHOLDS",scamp_sn_thresh,"-FWHM_THRESHOLDS",scamp_fwhm_thresh,"-CROSSID_RADIUS",scamp_crossid_rad,"-MATCH_RESOL",scamp_match_resol]
# 		command = ['/Users/mtaylor/local/bin/scamp','@'+scamp_list_file,'-c',SCAMP_config_dir+"/ctio_decam.scamp"]#,"-MERGEDOUTCAT_TYPE",scamp_cat_type_out,"-MERGEDOUTCAT_NAME",scamp_cat_file_out,"-MATCH","Y","-WRITE_XML","N","-SAVE_REFCATALOG","Y","-REFOUT_CATPATH",scamp_refcat_dir,"-CHECKPLOT_DEV","PSC","-CHECKPLOT_ANTIALIAS","Y","-CHECKPLOT_TYPE",scamp_check_type,"-CHECKPLOT_NAME",scamp_check_file,"-ASTREF_CATALOG",scamp_astref_cat,"-ASTREF_BAND",scamp_astref_band,"-ASTREFMAG_LIMITS",scamp_astref_maglim,"-DISTORT_DEGREES",scamp_distort_deg,"-PHOTCLIP_NSIGMA",scamp_photclip_nsig,"-SOLVE_ASTROM","Y","-POSITION_MAXERR",scamp_pos_err,"-PIXSCALE_MAXERR",scamp_scale_err,"-POSANGLE_MAXERR",scamp_angle_err,"-SN_THRESHOLDS",scamp_sn_thresh,"-CROSSID_RADIUS",scamp_crossid_rad,"-MATCH_RESOL",scamp_match_resol]
# 		command = ['/Users/mtaylor/local/bin/scamp','SE_output/sci_i_p1_d1_image_scamp.ldac','-c',SCAMP_config_dir+"/ctio_decam.scamp"]#,"-MERGEDOUTCAT_TYPE",scamp_cat_type_out,"-MERGEDOUTCAT_NAME",scamp_cat_file_out,"-MATCH","Y","-WRITE_XML","N","-SAVE_REFCATALOG","Y","-REFOUT_CATPATH",scamp_refcat_dir,"-CHECKPLOT_DEV","PSC","-CHECKPLOT_ANTIALIAS","Y","-CHECKPLOT_TYPE",scamp_check_type,"-CHECKPLOT_NAME",scamp_check_file,"-ASTREF_CATALOG",scamp_astref_cat,"-ASTREF_BAND",scamp_astref_band,"-ASTREFMAG_LIMITS",scamp_astref_maglim,"-DISTORT_DEGREES",scamp_distort_deg,"-PHOTCLIP_NSIGMA",scamp_photclip_nsig,"-SOLVE_ASTROM","Y","-POSITION_MAXERR",scamp_pos_err,"-PIXSCALE_MAXERR",scamp_scale_err,"-POSANGLE_MAXERR",scamp_angle_err,"-SN_THRESHOLDS",scamp_sn_thresh,"-CROSSID_RADIUS",scamp_crossid_rad,"-MATCH_RESOL",scamp_match_resol]
# 		for ii in range(len(command)):
# 			print command[ii], type(command[ii])
# 		exit()
		subprocess.call(command)
# 		exit()
# 						command = [SE_command,'temp_im_file.fits','-c',SE_config_file,'-CATALOG_NAME',SE_cat_file,'-PARAMETERS_NAME',SE_param_file,'-WEIGHT_IMAGE','temp_wt_file.fits','-XML_NAME',SE_xml_file,'-SATUR_LEVEL',str(SE_satur_level),'-MAG_ZEROPOINT',str(im_zp),'-SEEING_FWHM',str(im_fwhm),'-CHECKIMAGE_TYPE','BACKGROUND','-CHECKIMAGE_NAME',SE_check_file]
# 				subprocess.call(command)
