'''Use SWARP to create an image stack from images with
proper photometric and astrometric solutions.'''

def image_stack(filters,tiles,SWARP_config_dir,SWARP_out_dir,sci_im_dir,back_im_dir,stack_im_dir):
	import                 glob
	import                 os
	from astropy.io import fits
	import 		           subprocess
	
	swarp_command = "/Users/mtaylor/local/bin/swarp"

	#What kind of background subtraction will this be?
	#swarp = let swarp do the work
	#dynamic = background images have already been created; swarp will perform no subtraction
	back_type             = "swarp"
	swarp_backsize        = "256."

	swarp_config_file     = SWARP_config_dir+"ctio_decam.swarp"
	swarp_combine_type    = "CLIPPED"
	swarp_resampling_type = "LANCZOS2"
	swarp_resample_do     = "Y"
	swarp_resample_dir    = stack_im_dir+"/resample"
	if not os.path.exists(swarp_resample_dir): os.makedirs(swarp_resample_dir)

	for tile in tiles:
		for filter in filters:
			swarp_im_out = stack_im_dir+"/2014A-0610_FILTER%s_TILE%s_STACK.IMAGE.fits" % (filter,tile)
			swarp_wt_out = stack_im_dir+"/2014A-0610_FILTER%s_TILE%s_STACK.WEIGHT.fits" % (filter,tile)
	
			if back_type == "swarp":
				im_frameset = glob.glob("%s/sci_%s_*image.fits" % (sci_im_dir,filter))
				wt_frameset = [ii.replace("image","wtmap") for ii in im_frameset]

				#List of image/weight frames that will be assigned to frames with the bad chip removed.
				#These will the images that will actually be stacked.
				temp_im_list = ["temp_im_file%i.fits" % (im_number+1) for im_number in range(len(im_frameset))]
				temp_wt_list = [ii.replace("im","wt") for ii in temp_im_list]

# 				for jj in range(len(im_frameset)):
# 					print "Removing problem chip from %s and %s..." % (im_frameset[jj],wt_frameset[jj])
# 					#Remove chip S7 (extension 31) from calibration due to bad southern half of the chip
# 					temp_im_frame = fits.open(im_frameset[jj])
# 					temp_wt_frame = fits.open(wt_frameset[jj])
# 					im_frame = fits.HDUList()
# 					wt_frame = fits.HDUList()
# 					for ii in range(len(temp_im_frame)-1):
# 						if ii != 30:
# 							if ii < 30:
# 								im_frame.append(temp_im_frame[ii])
# 								wt_frame.append(temp_wt_frame[ii])
# 								im_frame[ii].header = temp_im_frame[ii].header
# 								wt_frame[ii].header = temp_wt_frame[ii].header
# 							else:
# 								im_frame.append(temp_im_frame[ii])
# 								wt_frame.append(temp_wt_frame[ii])
# 								im_frame[ii-1].header = temp_im_frame[ii].header
# 								wt_frame[ii-1].header = temp_wt_frame[ii].header
# 
# # 					temp_im_frame.close()
# # 					temp_wt_frame.close()
# 					im_frame.writeto(temp_im_list[jj],clobber=True)
# 					wt_frame.writeto(temp_wt_list[jj],clobber=True)
	
			if back_type == "dynamic":
				im_frameset = glob.glob("%s/sci_%s_*image.fits" % (back_im_dir,filter))
				wt_frameset = [ii.replace("image","wtmap") for ii in im_frameset]

		
			#Does the image list file exist? If not, create it. Otherwise add to it.
			swarp_list_file   = "%s/2014A-0610_FILTER%s_images.SWARP.list" % (SWARP_config_dir,filter)
			print "Creating SWARP image list: %s" % swarp_list_file
			temp_swarp_list_file = open(swarp_list_file,'w')
			for im_name in temp_im_list: print >> temp_swarp_list_file, im_name
			temp_swarp_list_file.close()

			command = [swarp_command,"@"+swarp_list_file,"-c","SWARP_config/default.swarp"]
			subprocess.call(command)



# 	swarp_head_file
# 	swarp_image_out
# 	swarp_weight_out
	
