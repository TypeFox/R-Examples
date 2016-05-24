resampleGridImagePoints <- function(pts, nx, rx, ry, fit.min.break=1, print.progress = FALSE){

	# CHECK IF ANY VALUES ARE NA
	if(is.na(pts[1,1])){
		if(print.progress) cat('Error:NA\n')
		return(list(pts=matrix(NA, nrow=rx*ry, ncol=2), error=rep(NA, rx*ry)))
	}

	# FIND SECOND GRID DIMENSION
	ny <- dim(pts)[1]/nx

	# MAX NUMBER OF TRIES WITH DIFFERENT INITIAL PARAMETERS
	max_fit_iter <- 15
	
	# LIST TO SAVE EACH FIT
	nlm_fits <- list()
	nlm_min <- rep(NA, max_fit_iter)

	for(i in 1:max_fit_iter){

		# SET INITIAL QUADRATIC CONSTANT PARAMETERS
		q_const <- rep((i-1)*(i-1)/max_fit_iter, 4)

		# SET INITIAL PARAMETERS
		p_init <- c(pts[1, 1], pts[1, 2], pts[nx, 1], pts[nx, 2], pts[nx*ny, 1], pts[nx*ny, 2], pts[(nx*ny)-nx+1, 1], pts[(nx*ny)-nx+1, 2], q_const)

		# FIND BEST FIT PARAMETERS
		#nlm_fits[[i]] <- nlm(imagePlaneGridTransformError, p=p_init, stepmax=100, iterlim=200, nx=nx, ny=ny, grid=pts)
		nlm_fits[[i]] <- nlminb(start=p_init, objective=imagePlaneGridTransformError, nx=nx, ny=ny, grid=pts)
		
		# SAVE MINIMUM
		nlm_min[i] <- nlm_fits[[i]]$objective
		
		# IF MINIMUM OF FIT IS LESS THAN fit.min.break, BREAK
		if(nlm_fits[[i]]$objective < fit.min.break) break
	}
	
	# SAVE FIT WITH LOWEST MINIMUM
	nlm_fit <- nlm_fits[[which.min(nlm_min)]]

	if(print.progress) cat('Mean fit error:', round(nlm_fit$objective, 4), 'px')

	#plot(pts)
	#points(imagePlaneGridTransform(nlm_fit$par, nx, ny), col='red')

	errors <- sqrt(rowSums((pts - imagePlaneGridTransform(nlm_fit$par, nx, ny))^2))
	if(print.progress) cat('; Max: ', round(max(errors), 4), ' px\n', sep='')

	pts <- imagePlaneGridTransform(nlm_fit$par, nx=rx, ny=ry)
	
	list(pts=pts, error=errors)
}