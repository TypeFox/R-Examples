estimateDistortion <- function(coor.2d, cal.nx, image.size){

	# GET NUMBER OF CORNERS IN OTHER DIMENSION
	cal.ny <- dim(coor.2d)[1] / cal.nx

	# SET PLANAR 3D COORDINATES (OBJECT COORDINATES)
	#coor_obj <- cbind(rep(0:(cal.nx-1), cal.ny), c(matrix(t(matrix((cal.ny-1):0, nrow=cal.ny, ncol=cal.nx)), nrow=1, byrow=F)), rep(0, cal.nx*cal.ny))
	coor_obj <- cbind(
		c(matrix(t(matrix(0:(cal.ny-1), nrow=cal.ny, ncol=cal.nx)), nrow=1, byrow=F)), 
		rep(0:(cal.nx-1), cal.ny), 
		rep(0, cal.nx*cal.ny)
	)

	# FILL ARRAY OF SAME DIMENSIONS AS 2D COORDINATES
	coor_obj_array <- array(coor_obj, dim = c(dim(coor_obj), dim(coor.2d)[3]))
	
	# REMOVE NA ASPECTS
	coor_obj_array <- coor_obj_array[, , !is.na(coor.2d[1,1,])]
	coor.2d <- coor.2d[, , !is.na(coor.2d[1,1,])]

	# TRY DIFFERENT STARTING PARAMETERS
	p_start <- list(
		rep(0, 5),
		c(0.1, NA, NA, NA, NA),
		c(0.1, 0.01, NA, NA, NA),
		c(0.14, 0.006, 0.002, NA, NA),
		c(0.1, 0.01, 1e-5, NA, NA),
		#c(0.1, 0.01, 1e-5, -1e-5, -0.1),
		c(0.01, -0.001, 1e-5, NA, NA),
		c(0.01, 0.1, -1e-5, NA, NA)
	)
	
	# SAVE WITH EACH TRY
	objectives <- rep(NA, length(p_start))
	par <- as.list(rep(NA, length(p_start)))
	
	# SAVE OBJECTIVE WITH NO DISTORTION
	objectives[1] <- distortionError(p=c(image.size[1]/2, image.size[2]/2, p_start[[1]]), 
		coor.img=coor.2d, coor.obj=coor_obj_array, image.size=image.size)
	par[[1]] <- c(image.size[1]/2, image.size[2]/2, p_start[[1]])

	# SKIP FIRST (NO DISTORTION CASE)
	for(i in 2:length(p_start)){

		# FIND OPTIMAL DISTORTION COEFFICIENTS, SKIP IF RETURNS ERROR
		nlm_fit <- tryCatch(
			expr={
				nlminb(start=c(image.size[1]/2, image.size[2]/2, p_start[[i]]), objective=distortionError, 
					coor.img=coor.2d, coor.obj=coor_obj_array, image.size=image.size)
			},
			error=function(cond) return(NULL),
			warning=function(cond) return(NULL)
		)

		if(is.null(nlm_fit)) next

		objectives[i] <- nlm_fit$objective
		par[[i]] <- nlm_fit$par
	}
	
	# GET PARAMETERS FROM RUN WITH LOWEST ERROR (INCLUDING NO DISTORTION CASE)
	dist_params <- par[[which.min(objectives)]]
	
	# ADD NAMES TO PARAMETERS
	names(dist_params) <- c('cx', 'cy', 'k1', 'k2', 'k3', 'p1', 'p2')
	
	dist_params
}