distortionError <- function(p, coor.img, coor.obj, image.size){

	p[is.na(p)] <- 0

	# APPLY DISTORTION PARAMETERS
	#coor_img_u <- undistort(coor.img, image.size=image.size, center=c(p[1], p[2]), 
	#	k=c(p[3], p[4], 0), p=c(p[5], p[6]))
	coor_img_u <- undistort(coor.img, image.size=image.size, center=c(p[1], p[2]), 
		k=c(p[3], p[4], p[5]), p=c(p[6], p[7]))
	
	# FIND HOMOGRAPHY ERROR
	find_hom <- findHomography(coor_img_u, coor.obj)

	mean(find_hom$error)
}
