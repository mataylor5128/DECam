'''Derive the photometric calibration for a set of
standard star frames taken with DECam'''

#----------------------------------------------------------#
#Compute Zeropoint, Airmass term, and Colour term Estimates#
#----------------------------------------------------------#
def photo_calib(calib_im_dir,calib_out_dir,sci_im_dir,SE_command,SE_config_dir,SE_output_dir,filters):
	import numpy              as     np
	import matplotlib.pyplot  as     plt
	from astropy.io           import fits, ascii
	from astropy              import units as u
	from astropy.time         import Time
	from astropy.coordinates  import SkyCoord, EarthLocation, AltAz
	from astropy.stats        import sigma_clip
	import                           os
	import                           subprocess
	import                           sklearn
	from sklearn.linear_model import LinearRegression
	import time               as     tm

	#Open new calibration output file
	calib_file_out = "%s/scabs_2014A-0610_calibration-%s.txt" % (calib_out_dir,Time(tm.strftime("%Y-%m-%d %I:%M:%S")).mjd)

	if not os.path.exists(calib_file_out) or overwrite == "T":
		calib_data_out = open(calib_file_out,'w')
		print >> calib_data_out, "#filter	MJD-OBS		zp	zp_err	k	k_err	m0	m0_err	a	a_err"
	else:
		calib_data_out = open(calib_file_out,'a')	

	#Loop through each requested filter
	for filter in filters:
		#Set up SE input/output
		SE_config_file = SE_config_dir+"/ctio_decam_zp.sex"
		SE_param_file  = SE_config_dir+"/ctio_decam_zp.param"
		SE_filt_name   = SE_config_dir+"/gauss_3.5_7x7.conv"

		if filter == 'z':
			im_frameset = ["%s/std-lse_44-%s_2s_image_1.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_2s_image_3.fits" % (calib_im_dir,filter),# "%s/std-lse_44-%s_2s_image_3.fits" % (calib_im_dir,filter), #"%s/std-lse_44-%s_2s_image_4.fits" % (calib_im_dir,filter),
						"%s/std-lse_44-%s_10s_image_1.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_10s_image_3.fits" % (calib_im_dir,filter),# "%s/std-lse_44-%s_10s_image_3.fits" % (calib_im_dir,filter), #"%s/std-lse_44-%s_10s_image_4.fits" % (calib_im_dir,filter),
						"%s/std-lse_44-%s_30s_image_1.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_30s_image_3.fits" % (calib_im_dir,filter)]#, "%s/std-lse_44-%s_30s_image_3.fits" % (calib_im_dir,filter)]#, "%s/std-lse_44-%s_30s_image_4.fits" % (calib_im_dir,filter)]
		else:	
			im_frameset = ["%s/std-lse_44-%s_2s_image_1.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_2s_image_2.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_2s_image_3.fits" % (calib_im_dir,filter), #"%s/std-lse_44-%s_2s_image_4.fits" % (calib_im_dir,filter),
						"%s/std-lse_44-%s_10s_image_1.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_10s_image_2.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_10s_image_3.fits" % (calib_im_dir,filter), #"%s/std-lse_44-%s_10s_image_4.fits" % (calib_im_dir,filter),
						"%s/std-lse_44-%s_30s_image_1.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_30s_image_2.fits" % (calib_im_dir,filter), "%s/std-lse_44-%s_30s_image_3.fits" % (calib_im_dir,filter)]#, "%s/std-lse_44-%s_30s_image_4.fits" % (calib_im_dir,filter)]
		wt_frameset = [ii.replace("image","wtmap") for ii in im_frameset]

		#Initialize the output data frame array
		zero_data = []

		for ii in range(len(im_frameset)):
			print "Processing frame: %s" % im_frameset[ii]
			#Get basic image info
		# 	im_name = "%s/%s" % im_frameset#(calib_im_dir, im_frameset[ii])
		# 	wt_name = "%s/%s" % im_frameset#(calib_im_dir, wt_frameset[ii])
			im_name = im_frameset[ii]#"%s/%s" % im_frameset#(calib_im_dir, im_frameset[ii])
			wt_name = wt_frameset[ii]#"%s/%s" % im_frameset#(calib_im_dir, wt_frameset[ii])

			SE_cat_file = im_name.replace(calib_im_dir,SE_output_dir).replace(".fits","_zp.ldac")
			SE_xml_file = im_name.replace(calib_im_dir,SE_output_dir).replace(".fits","_zp.xml")
			SE_check_file = im_name.replace(calib_im_dir,SE_output_dir).replace(".fits",".BACK.fits")

			#Open the image and weightmap
			temp_im_frame = fits.open(im_name)
			temp_wt_frame = fits.open(wt_name)

			#Get some basic header information
			im_h = temp_im_frame[0].header
			im_filter = im_h['FILTER'][0]
			im_exptime = im_h['EXPTIME']
			im_zp = im_h['MAGZERO']
			im_coord = SkyCoord(im_h['RA'],im_h['DEC'], unit=(u.hourangle, u.deg))
			im_ra = im_coord.ra.deg
			im_dec = im_coord.dec.deg
			im_mjd = im_h['MJD-OBS']

			#Calculate airmass of pointing
			ctio = EarthLocation(lat=im_h['OBS-LAT']*u.deg,lon=-1*im_h['OBS-LONG']*u.deg,height=im_h['OBS-ELEV']*u.m)
			time = Time(im_h['DATE-OBS'].replace('T',' '))
			im_altaz = im_coord.transform_to(AltAz(obstime=time,location=ctio,))
			im_airmass = im_altaz.secz

			#Has the photometry for the zeropoint estimation already been extracted?
			#If not, proceed to extract photometry...
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
				command = [SE_command,'temp_im_file.fits','-c',SE_config_file,'-CATALOG_NAME',SE_cat_file,'-PARAMETERS_NAME',SE_param_file,'-WEIGHT_IMAGE','temp_wt_file.fits','-XML_NAME',SE_xml_file,'-SATUR_LEVEL',str(SE_satur_level),'-MAG_ZEROPOINT',str(im_zp),'-SEEING_FWHM',str(im_fwhm),'-CHECKIMAGE_TYPE','BACKGROUND','-CHECKIMAGE_NAME',SE_check_file]
				subprocess.call(command)

			#Free up memory
			temp_im_frame.close()
			temp_wt_frame.close()

			#Open the SE catalogue file and extract source coordinates from each chip into
			#a single catalogue.
			SE_data_hdu = fits.open(SE_cat_file)
			ii = 1
			while ii < len(SE_data_hdu):
				if ii == 1:
					SE_data = SE_data_hdu[ii+1].data
				else:
					SE_data = np.concatenate([SE_data,SE_data_hdu[ii+1].data])
				ii += 2

			#Eliminate sources with FLAGS > 0
			SE_data = SE_data[(SE_data['FLAGS'] == 0)]
			
			SE_source_coords = SkyCoord(SE_data['ALPHA_J2000'],SE_data['DELTA_J2000'],unit=(u.deg,u.deg))

			#Open the catalogue of standard stars and compile the list of successful matches
			print "Matching standard stars..."
			standard_type = "southern" # "sdss"
			if standard_type == "southern":
				std_cat_file = "/Users/mtaylor/Projects/SCABS/CN2014A-92/data/standard_stars/Southern_ugriz_trimmed/LSE_44.dat.trimmed"

			std_cat_data = ascii.read(std_cat_file)

			std_cat_coords = SkyCoord(std_cat_data['col2'], std_cat_data['col3'], unit=(u.hourangle,u.deg))

			idx, d2d,d3d = std_cat_coords.match_to_catalog_sky(SE_source_coords)

			#Limit matches to only those with separations < 0.5 arcsec
			mask = (d2d < 0.25*u.arcsec)
			d2d = d2d[mask]
			idx = idx[mask]
			std_cat_match = std_cat_data[mask]
			SE_data_match = [SE_data[ii] for ii in idx]
			SE_data_match = np.array(SE_data_match)

			std_match_coords = SkyCoord(std_cat_match['col2'],std_cat_match['col3'],unit=(u.hourangle,u.deg))
			SE_match_coords  = SkyCoord(SE_data_match['ALPHA_J2000'],SE_data_match['DELTA_J2000'],unit=(u.deg,u.deg))

			if filter == 'u':
				mag_std = std_cat_match['col4']
			if filter == 'g':
				mag_std = std_cat_match['col7']
			if filter == 'r':
				mag_std = std_cat_match['col10']
			if filter == 'i':
				mag_std = std_cat_match['col13']
			if filter == 'z':
				mag_std = std_cat_match['col16']
		# 	mag_SE = SE_data_match['MAG_APER'][:,4]+2.5*np.log10(im_exptime)-im_zp
			mag_SE  = SE_data_match['MAG_AUTO']+2.5*np.log10(im_exptime)-im_zp
			flux_rad = SE_data_match['FLUX_RADIUS']
			airmass = [im_airmass]*len(mag_std)
			exptime = [im_exptime]*len(mag_std)

			if filter == "u":
				std_col = std_cat_match['col4']-std_cat_match['col13']
				color_label = "u' - i'"
			if filter == "g":
				std_col = std_cat_match['col7']-std_cat_match['col10']
				color_label = "g' - r'"
			if filter == "r":
				std_col = std_cat_match['col7']-std_cat_match['col10']
				color_label = "g' - r'"
			if filter == "i":
				std_col = std_cat_match['col13']-std_cat_match['col16']
				color_label = "i' - z'"
			if filter == "z":
				std_col = std_cat_match['col13']-std_cat_match['col16']
				color_label = "i' - z'"

			zero_data.append({"mag_std": np.array(mag_std), "mag_SE": np.array(mag_SE),
								"delta_mag": np.array(mag_std) - np.array(mag_SE), 
								"airmass": np.array(airmass), 
								"std_color": np.array(std_col),
								"exptime": np.array(exptime),
								"flux_rad": np.array(flux_rad)})

		#sigma-clip each set of zeropoint data
		for ii in range(len(zero_data)):
			mask = sigma_clip(zero_data[ii]['delta_mag'],sig=1.).mask
			zero_data[ii]["mag_std"] = zero_data[ii]["mag_std"][~mask]
			zero_data[ii]["mag_SE"] = zero_data[ii]["mag_SE"][~mask]
			zero_data[ii]["delta_mag"] = zero_data[ii]["delta_mag"][~mask]
			zero_data[ii]["airmass"] = zero_data[ii]["airmass"][~mask]
			zero_data[ii]["std_color"] = zero_data[ii]["std_color"][~mask]
			zero_data[ii]["flux_rad"] = zero_data[ii]["flux_rad"][~mask]

		med_dmag = []
		med_airmass = []
		med_color = []
		for ii in range(len(zero_data)):
			med_dmag.append(np.median(zero_data[ii]['delta_mag']))
			med_airmass.append(np.median(zero_data[ii]['airmass']))
			med_color.append(np.median(zero_data[ii]['std_color']))

		for ii in range(len(zero_data)):
			if ii == 0:
				mag_std = zero_data[ii]['mag_std']
				mag_SE = zero_data[ii]['mag_SE']
				dmag = zero_data[ii]['delta_mag']
				airmass = zero_data[ii]['airmass']
				color = zero_data[ii]['std_color']
				flux_rad = zero_data[ii]['flux_rad']
			else:
				mag_std = np.concatenate([mag_std,zero_data[ii]['mag_std']])
				mag_SE = np.concatenate([mag_SE,zero_data[ii]['mag_SE']])
				dmag = np.concatenate([dmag,zero_data[ii]['delta_mag']])	
				airmass = np.concatenate([airmass,zero_data[ii]['airmass']])
				color = np.concatenate([color,zero_data[ii]['std_color']])
				flux_rad = np.concatenate([flux_rad,zero_data[ii]['flux_rad']])

		#Derive the airmass term and zeropoint using a bootstrap approach
		#Apply linear regression model to the dmag v. airmass plot
		#Adopt slope as airmass term and intercept as zeropoint

		#Initialize output arrays for k and zp
		k_out = []
		zp_out = []

		#Plot airmass vs. delta_mag for the standard star observations
		plt.figure()
		for ii in range(len(zero_data)):
			if zero_data[ii]['exptime'][0] == 2:
				plt.plot(zero_data[ii]['airmass'],zero_data[ii]['delta_mag'],'bo',mec=None,alpha=0.5,label='2s')
			if zero_data[ii]['exptime'][0] == 10:
				plt.plot(zero_data[ii]['airmass'],zero_data[ii]['delta_mag'],'ro',mec=None,alpha=0.5,label='10s')
			if zero_data[ii]['exptime'][0] == 30:
				plt.plot(zero_data[ii]['airmass'],zero_data[ii]['delta_mag'],'go',mec=None,alpha=0.5,label='30s')

		#Set up a linear regression model
		xx = np.linspace(0,2.5,1e3)
		regr = LinearRegression(fit_intercept=True)
		indices = np.arange(len(dmag))

		#Resample the true data with replacement n_trials times, and apply
		#a linear regression each time. Record the best fit parameters for each
		#iteration.
		N = 0
		n_trials = 1000
		while N < n_trials:
			#resample indices and select corresponding delta_mag and airmasses
			temp_indices = np.random.choice(indices,size=len(indices),replace=True)
			temp_dmag = np.array([dmag[ind] for ind in temp_indices]).reshape((-1,1))
			temp_airmass = np.array([airmass[ind] for ind in temp_indices]).reshape((-1,1))

			#apply linear regression and record the parameters
			regr.fit(temp_airmass,temp_dmag)
			temp_k = regr.coef_[0]
			temp_zp = regr.intercept_[0]

			k_out.append(temp_k[0])
			zp_out.append(temp_zp)
 
			#lightly overplot the fit for this iteration
			yy = temp_k*xx+temp_zp
			plt.plot(xx,yy,'k--',alpha=0.1,lw=0.1)

			#iterate
			N += 1

		#Adopt the mean of the the outputs as the final values,
		#and 1sigma as the errors.
		k_err = np.std(k_out)
		k_out = np.mean(k_out)
		zp_err = np.std(zp_out)
		zp_out = np.mean(zp_out)
		print "Airmass term = %.2f +/- %.2f" % (k_out,k_err)
		print "Zeropoint = %.2f +/- %.2f" % (zp_out,zp_err)

		#Overplot the adopted airmass correction function
		yy = k_out*xx+zp_out
		plt.plot(xx,yy,'k--',lw=3.)

		#Finalize plot
		plt.annotate(r"$\Delta %s'  = (%.2f\pm%.2f)\times X + (%.2f\pm%.2f)$" % (filter,k_out,k_err,zp_out,zp_err),xy=(0.15,0.9),xycoords='axes fraction',alpha=0.5,fontsize=20)
		plt.annotate(r'2s',xy=(0.05,0.15),xycoords='axes fraction',color='b',alpha=0.5,fontsize=20)
		plt.annotate(r'10s',xy=(0.05,0.1),xycoords='axes fraction',color='r',alpha=0.5,fontsize=20)
		plt.annotate(r'30s',xy=(0.05,0.05),xycoords='axes fraction',color='g',alpha=0.5,fontsize=20)
		plt.xlabel(r'AIRMASS',fontsize=24)
		plt.ylabel(r"$%s'_{\rm std} - %s'_{\rm inst}$ [mag]" % (filter,filter),fontsize=24)
		plt.xlim(0,2.2)
		if filter != 'u':
			plt.ylim(med_dmag[0]-0.5,med_dmag[0]+0.5)
		else:
			plt.ylim(med_dmag[0]-0.5,med_dmag[0]+1.0)
		plt.savefig("%s/dmag_v_airmass_mjd%i_%s.pdf" % (calib_out_dir,int(im_mjd),filter), bbox_inches="tight")

		#Correct m_inst for atmospheric extinction, and isolate 
		#the best delta_mag matches to derive the colour term
		#by clipping to within +/- 1sigma of the mean
		dmag -= zp_out + k_out*airmass

		#Plot the total sample
		plt.figure()
		plt.plot(color,dmag,'ro')

		#Apply sigma-clipping
		mask = sigma_clip(dmag,sig=1.).mask
		dmag = dmag[~mask]
		color = color[~mask]
		airmass = airmass[~mask]

		#As with the airmass corrections, use a bootstrapping
		#approach to derive the colour correction

		#Initialize output arrays for a and m_0
		a_out = []
		m0_out = []

		#Set up a linear regression model
		xx = np.linspace(np.min(color)-0.05,np.max(color)+0.05,1e3)
		regr = LinearRegression(fit_intercept=True)
		indices = np.arange(len(dmag))

		#Resample the true data with replacement n_trials times, and apply
		#a linear regression each time. Record the best fit parameters for each
		#iteration.
		N = 0
		n_trials = 1000
		while N < n_trials:
			#resample indices and select corresponding delta_mag and colors
			temp_indices = np.random.choice(indices,size=len(indices),replace=True)
			temp_dmag = np.array([dmag[ind] for ind in temp_indices]).reshape((-1,1))
			temp_color = np.array([color[ind] for ind in temp_indices]).reshape((-1,1))

			#apply linear regression and record the parameters
			regr.fit(temp_color,temp_dmag)
			temp_a = regr.coef_[0]
			temp_m0 = regr.intercept_[0]

			a_out.append(temp_a[0])
			m0_out.append(temp_m0)
 
			#lightly overplot the fit for this iteration
			yy = temp_a*xx+temp_m0
			plt.plot(xx,yy,'k--',alpha=0.1,lw=0.1)

			#iterate
			N += 1

		#Adopt the mean of the the outputs as the final values,
		#and 1sigma as the errors.
		a_err = np.std(a_out)
		a_out = np.mean(a_out)
		m0_err = np.std(m0_out)
		m0_out = np.mean(m0_out)
		print "Colour term = %.2f +/- %.2f" % (a_out,a_err)
		print "Zeropoint = %.2f +/- %.2f" % (m0_out,m0_err)

		#Overplot the adopted airmass correction function
		yy = a_out*xx+m0_out
		plt.plot(xx,yy,'k--',lw=3.)

		#Finalize plot
		plt.plot(color,dmag,'bo')
		plt.plot(xx,yy,'k--',lw=2)
		plt.annotate(r"$\Delta %s'  = (%.2f\pm%.2f)\times (%s) + (%.2f\pm%.2f)$" % (filter,a_out,a_err,color_label,m0_out,m0_err),xy=(0.15,0.9),xycoords='axes fraction',alpha=0.5,fontsize=20)
		if filter == 'u':
			plt.ylim(-0.5,0.5)
		plt.xlabel(r"$(%s)$ [mag]" % color_label,fontsize=24)
		plt.ylabel(r"$%s'_{\rm std} - %s'_{\rm inst}$ [mag]" % (filter,filter),fontsize=24)
		plt.grid()
		plt.savefig("%s/dmag_v_color_mjd%i_%s.pdf" % (calib_out_dir,int(im_mjd),filter), bbox_inches="tight")

		#Correct magnitudes for the colour term/zeropoint
		dmag -= m0_out + a_out*color

		#Plot the final calibrated magnitudes
		plt.figure()
		plt.subplot(211)
		plt.plot(airmass,dmag,'bo')
		plt.xlabel(r'AIRMASS',fontsize=24)
		plt.ylabel(r"$%s'_{\rm std} - %s'_{\rm inst}$ [mag]" % (filter,filter),fontsize=24)
		plt.grid()
		plt.subplot(212)
		plt.plot(color,dmag,'bo')
		plt.xlabel(r"$(%s)$ [mag]" % color_label,fontsize=24)
		plt.ylabel(r"$%s'_{\rm std} - %s'_{\rm inst}$ [mag]" % (filter,filter),fontsize=24)
		plt.grid()
		plt.savefig("%s/final-calib_mjd%i_%s.pdf" % (calib_out_dir,int(im_mjd),filter), bbox_inches="tight")

		#Write results to file
		print >> calib_data_out, "%s	%f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f	%.2f" % (filter,im_mjd,zp_out,zp_err,k_out,k_err,m0_out,m0_err,a_out,a_err)

	calib_data_out.close()
	subprocess.call(['rm','temp_im_file.fits'])
	subprocess.call(['rm','temp_wt_file.fits'])

#Saturation = 65000 ADU

#SE for SCAMP input

#SCAMP for astrometric solution

#Background subtraction

#SWARP for final image stack