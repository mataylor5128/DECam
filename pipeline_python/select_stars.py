'''Algorithm to isolate stellar sources
from other sources on a CCD chip using SE's
MAG_AUTO and FLUX_RADIUS parameters'''
import numpy as np
import numpy.ma as ma
from astropy.stats import sigma_clip
import matplotlib.pyplot as plt

def find_mag_lim(med_flux_rads,mean_mags):
	ii = 1
	#Find brightest mag_lim
	#Check flux_rad for points immediately before/after mag_chk, are they similar to flux_chk?
	while ii < len(med_flux_rads):
		flux_chk = med_flux_rads[ii]
		flux_low = med_flux_rads[ii-1]
		flux_high = med_flux_rads[ii+1]
		dflux_low = np.abs(flux_chk - flux_low)
		dflux_high = np.abs(flux_chk - flux_high)
		if dflux_low <= 0.05 and dflux_high <= 0.05:
			mag_lim_low = mean_mags[list(med_flux_rads).index(flux_chk)]
			ii = len(med_flux_rads)
		else:
			ii += 1
	ii = len(med_flux_rads) - 2
	#Find faintest mag_lim
	#Same as bright, but working backwards
	while ii > 0:
		flux_chk = med_flux_rads[ii]
		flux_low = med_flux_rads[ii-1]
		flux_high = med_flux_rads[ii+1]
		dflux_low = np.abs(flux_chk - flux_low)
		dflux_high = np.abs(flux_chk - flux_high)
		if dflux_low <= 0.25 and dflux_high <= 0.25:
			mag_lim_high = mean_mags[list(med_flux_rads).index(flux_chk)]
			ii = 0
		else:
			ii -= 1

	return mag_lim_low, mag_lim_high
	

def select_stars(flux_rad, mag_auto):	
	'''Given a set of SE parameters flux_radius and mag_auto, isolate the locus
	of stellar sources, and pick appropriate, non-saturated stars.'''

	plotting = False

	#Break the magnitudes and fluxes in bins and calculate statistics in each bin
	mag_chk = np.min(mag_auto)
	dmag = 0.25
	dmag_max = 2.5

	mean_mags = []
	mean_flux_rads = []
	median_flux_rads = []
	stddev_flux_rads = []
	while mag_chk < np.max(mag_auto):
		#Is dmag > dmag_max and no sources in sample? If so, end the loop
		if dmag > dmag_max and len(flux_rad[np.logical_and(mag_auto >= mag_chk, mag_auto < mag_chk + dmag)]) < 2:
			mag_chk = np.max(mag_auto)
		else:
			#Are there > 1 sources in the sample?
			if len(flux_rad[np.logical_and(mag_auto >= mag_chk, mag_auto < mag_chk + dmag)]) < 2:
				dmag += 0.05
			else:
				temp_flux_rad = flux_rad[np.logical_and(mag_auto >= mag_chk, mag_auto < mag_chk + dmag)]
				temp_mag_auto = mag_auto[np.logical_and(mag_auto >= mag_chk, mag_auto < mag_chk + dmag)]

				mean_mags.append(np.mean(temp_mag_auto))
				mean_flux_rads.append(np.mean(temp_flux_rad))
				median_flux_rads.append(np.median(temp_flux_rad))

				dmag = 0.25			
				mag_chk += dmag

	mean_mags = np.array(mean_mags)	
	mean_flux_rads = np.array(mean_flux_rads)
	median_flux_rads = np.array(median_flux_rads)	

	if plotting == True:
		fig = plt.figure()
		ax = plt.subplot()
		plt.plot(mean_mags,mean_flux_rads,'o')
		plt.plot(mean_mags,median_flux_rads,'^')
		plt.xlabel(r"MAG_AUTO")
		plt.ylabel(r"FLUX_RADIUS")
		plt.savefig("meanmag_v_fluxrad.pdf",bbox_inches="tight")
		plt.show()
		exit()

	#Determine bright and faint magnitude limits using the statistics in each bin
	mag_lim_low, mag_lim_high = find_mag_lim(median_flux_rads,mean_mags)

	#Break the magnitude limits into 'nbin' bins
	nbins = 5
	dmag = (mag_lim_high - mag_lim_low)/nbins
	mag_chk = mag_lim_low
	array_out = [False] * len(flux_rad)
	array_out = np.array(array_out)
	while mag_chk < mag_lim_high:

		#---------------------------------#
		#Select stars within magnitude bin#
		#---------------------------------#
		#Mask sources within the magnitude bin.
		temp_mag_mask = ma.masked_outside(mag_auto, mag_chk, mag_chk+dmag)
		gv_inds = np.where(temp_mag_mask.mask == False)
		
		#Get flux_radius values corresponding to gv_inds
		temp_flux_rads = [flux_rad[ii] for ii in gv_inds]
		
		#Sigma-clip outliers from magnitude limited subsample, record the indices from gv_inds,
		#and recover clipped, magnitude limited flux_radii
		clip_flux_rad = sigma_clip(temp_flux_rads,sig=2.3)
		gv_clipped_inds = [gv_inds[0][ii] for ii in np.where(clip_flux_rad.mask[0] == False)[0]]
		temp_flux_rads = [flux_rad[ii] for ii in gv_clipped_inds]

		#Calculate dispersion of clipped, magnitude limited flux_radii and record indices of 
		#flux_radius corresponding to within +/- 1-sigma of the median flux
		temp_stddev = np.std(temp_flux_rads)
		cut_factor = 1.*temp_stddev
		temp_flux_rad_mask = np.logical_and(temp_flux_rads >= np.median(temp_flux_rads)-cut_factor, temp_flux_rads <= np.median(temp_flux_rads)+cut_factor)
		gv_star_inds = np.array(gv_clipped_inds)[temp_flux_rad_mask]
		
		#If this is the 1st magnitude bin, initialize output index array, otherwise
		#concatenate to previous output index array. Increase mag_chk to the next bin.
		if mag_chk == mag_lim_low:
			gv_star_inds_out = gv_star_inds
		else:
			gv_star_inds_out = np.concatenate([gv_star_inds_out,gv_star_inds])
		mag_chk += dmag

		#Check for duplicates
		dupes = set([x for x in gv_star_inds_out if list(gv_star_inds_out).count(x) > 1])
		if len(dupes) > 0:
			print dupes
			print "Duplicate stellar sources found! Exiting..."
			exit()

	return gv_star_inds_out