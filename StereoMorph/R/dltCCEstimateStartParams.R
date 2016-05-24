dltCCEstimateStartParams <- function(coor.2d, num.grid, nx, ny, grid.size, p.fix, use.param.min=3){

	nrow_start_param <- 576

	set_seed <- FALSE
	
	if(!is.null(p.fix) && length(p.fix) >= use.param.min*6){

		# GET 3D COORDINATES BASED ON OPTIMIZED PARAMETERS
		coor_3d_coeff <- transformPlanarCalibrationCoordinates(tpar=c(matrix(c(rep(0, 6), p.fix), nrow=6)), 
			nx=nx, ny=ny, sx=grid.size)

		# GET 2D INPUT COORDINATES
		coor_2d_coeff <- apply(coor.2d[, , 1:(num.grid-1), ], c(2, 4), matrix, byrow=FALSE)
		
		# GET DLT CALIBRATION COEFFICIENTS FROM OPTIMIZED 3D COORDINATE SUBSET (ALL GRID POINTS - SO RMSE CAN DIFFER FROM NLM MINIMUM)
		dlt_coefficients_t <- dltCoefficients(coor.3d=coor_3d_coeff, coor.2d=coor_2d_coeff)
		
		# USE DLT COEFFICIENTS TO RECONSTRUCT ALL 2D COORDINATES
		dlt_reconstruct <- dltReconstruct(cal.coeff=dlt_coefficients_t$cal.coeff, coor.2d=apply(coor.2d[, , 1:num.grid, ], c(2, 4), matrix, byrow=FALSE))

		# COPY 3D COORDINATES TO ARRAY
		grids_3d <- array(NA, dim=c(nx*ny, 3, num.grid))
		for(i in 1:dim(grids_3d)[3]) grids_3d[, , i] <- dlt_reconstruct$coor.3d[((i-1)*nx*ny+1):(i*nx*ny), ]
		
		# FIND INVERSE TRANSFORMATIONS TO GET 
		t_param <- inverseGridTransform(grids_3d)
		
		# STARTING PARAMETER ARRAY
		start_param_array <- array(NA, dim=c(nrow_start_param, num.grid-1, 6))
		for(i in 1:dim(start_param_array)[2]) start_param_array[, i, ] <- matrix(t_param[i+1, ], nrow=nrow_start_param, ncol=6, byrow=TRUE)

		# APPLY RANDOM SCALING
		if(set_seed) set.seed(42);
		scaling <- array(sample(seq(0.5, 1.5, length=dim(start_param_array)[1]*dim(start_param_array)[2]*3)), dim=c(dim(start_param_array)[1], dim(start_param_array)[2], 3))
		start_param_array[, , 4:6] <- array(start_param_array[, , 4:6], dim=dim(scaling))*scaling

		# APPLY RANDOM SCALING
		if(set_seed) set.seed(42);
		scaling <- array(sample(seq(0, 2, length=dim(start_param_array)[1]*dim(start_param_array)[2]*3)), dim=c(dim(start_param_array)[1], dim(start_param_array)[2], 3))
		start_param_array[, , 1:3] <- array(start_param_array[, , 1:3], dim=dim(scaling))*scaling

		max_disp_rwu <- max(abs(t_param[, 4:6]))

	}else{

		# LIMIT TO NUMBER OF GRIDS IN CONSIDERATION
		coor.2d <- coor.2d[, , 1:num.grid, ]

		# GET APPROXIMATE SCALING FROM PIXELS TO REAL-WORLD UNITS
		px2rwu <- 0
		for(i in 1:dim(coor.2d)[4]){
			for(j in 1:dim(coor.2d)[3]){
				dimm <- (grid.size*(nx-1)) / sqrt(sum((coor.2d[1, , j, i] - coor.2d[nx, , j, i])^2))
				if(is.na(dimm)) next
				px2rwu <- max(px2rwu, (grid.size*(nx-1)) / sqrt(sum((coor.2d[1, , j, i] - coor.2d[nx, , j, i])^2)))
			}
		}

		# DISTANCE OF CENTROIDS OF ALL GRID ASPECTS FROM FIRST ASPECT CENTROID
		centroid_dist <- matrix(NA, dim(coor.2d)[3]-1, dim(coor.2d)[4])
		for(i in 1:dim(coor.2d)[4]){

			# CENTROID OF FIRST ASPECT
			centroid_first <- colMeans(coor.2d[, , 1, i])
		
			# CENTROID-CENTROID DISTANCE
			for(j in 2:dim(coor.2d)[3]) centroid_dist[j-1, i] <- sqrt(sum((colMeans(coor.2d[, , j, i]) -  centroid_first)^2))
		}
		
		# REPLACE NA VALUES WITH MEAN
		centroid_dist[is.na(centroid_dist)] <- mean(centroid_dist, na.rm=TRUE)

		# FIND MAXIMUM CENTROID DISPLACEMENT IN REAL-WORLD UNITS
		max_disp_rwu <- max(centroid_dist, na.rm=TRUE)*px2rwu

		# FIND MAX SIDE DIMENSIONS IN STANDARD SQUARE LENGTH (RELATIVE MEASURE OF DISTANCE FROM CAMERA)
		square_side_max <- matrix(NA, dim(coor.2d)[3], dim(coor.2d)[4])
		for(i in 1:dim(coor.2d)[4]){
			for(j in 1:dim(coor.2d)[3]){
			
				# SIDE LENGTHS
				side_lengths <- c(
					sqrt(sum((coor.2d[1, , j, i] - coor.2d[nx, , j, i])^2, na.rm=TRUE)),
					sqrt(sum((coor.2d[nx*ny-nx+1, , j, i] - coor.2d[nx*ny, , j, i])^2, na.rm=TRUE)),
					sqrt(sum((coor.2d[nx, , j, i] - coor.2d[nx*ny, , j, i])^2, na.rm=TRUE)),
					sqrt(sum((coor.2d[1, , j, i] - coor.2d[nx*ny-nx+1, , j, i])^2, na.rm=TRUE)))
				
				# FIND MAX AND SCALE
				if(which.max(side_lengths) %in% c(1, 2)){
					square_side_max[j, i] <- max(side_lengths)
				}else{
					square_side_max[j, i] <- max(side_lengths)*(nx / ny)
				}
			}
		}

		# FIND DIFFERENCE FROM FIRST ASPECT
		square_side_dist <- square_side_max[2:nrow(square_side_max), ] - matrix(square_side_max[1, ], nrow(square_side_max)-1, ncol(square_side_max), byrow=TRUE)

		# SCALE TO MAX OF CENTROID DISTANCES
		square_side_dist <- abs(square_side_dist) * (max(centroid_dist, na.rm=TRUE) / max(abs(square_side_dist), na.rm=TRUE))

		# CONVERT TO REAL-WORLD UNITS AND COMBINE INTO ONE MATRIX
		dist_param <- cbind(centroid_dist * px2rwu, square_side_dist * px2rwu)

		# STARTING PARAMETER ARRAY
		start_param_array <- array(NA, dim=c(nrow_start_param, dim(coor.2d)[3]-1, 6))

		# ADD IN TRANSLATION STARTING PARAMETERS
		for(i in 1:dim(start_param_array)[2]){
			start_param_array[, i, 4:6] <- c(dist_param[i, 1:4], 0, dist_param[i, 4:1], dist_param[i, 3:2], dist_param[i, 1], 0, dist_param[i, 4], dist_param[i, 1:2], 0, dist_param[i, 4])
			start_param_array[, i, 4:6] <- start_param_array[, i, 4:6]*c(-1,1,1,-1)
		}

		# ADD IN ROTATION STARTING PARAMETERS
		start_param_array[, , 1:3] <- c(0.0216, 0.0113, 0.0124, 0.02, 0.0348, 0.0323, 0.0295, 0.02143, 0.034)
		start_param_array[, , 1:3] <- start_param_array[, , 1:3]*c(-1,1,1,-1)

		# APPLY RANDOM SCALING
		if(set_seed) set.seed(42);
		scaling <- array(sample(seq(0.7, 1.4, length=dim(start_param_array)[1]*dim(start_param_array)[2]*3)), dim=c(dim(start_param_array)[1], dim(start_param_array)[2], 3))
		start_param_array[, , 4:6] <- array(start_param_array[, , 4:6], dim=dim(scaling))*scaling

		# APPLY RANDOM SCALING
		if(set_seed) set.seed(42);
		scaling <- array(sample(seq(0.7, 2, length=dim(start_param_array)[1]*dim(start_param_array)[2]*3)), dim=c(dim(start_param_array)[1], dim(start_param_array)[2], 3))
		start_param_array[, , 1:3] <- array(start_param_array[, , 1:3], dim=dim(scaling))*scaling
	}

	#cat('\n')
	#print(start_param_array[, 1, 1:3])

	list(start_param_array = start_param_array, max_disp_rwu = max_disp_rwu)
}