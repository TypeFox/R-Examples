dltReconstruct <- function(cal.coeff, coor.2d, min.views = 2){
	# Translated and modified from Matlab function dlt_reconstruct() written by Ty Hedrick

	# IF INPUT IS LIST, CONVERT TO MATRIX
	if(is.list(coor.2d)) coor.2d <- landmarkListToMatrix(coor.2d)
	
	# IF INPUT IS ARRAY, CONVERT TO MATRIX
	if(length(dim(coor.2d)) == 3) coor.2d <- t(apply(coor.2d, 1, matrix))		

	# MATRIX FOR POINT RECONSTRUCTION
	coor.3d <- matrix(NA, nrow=nrow(coor.2d), ncol=3, dimnames=list(rownames(coor.2d), NULL))

	# INTERNAL ROOT MEAN SQUARE ERROR
	rmse <- rep(NA, nrow(coor.2d))

	for(i in 1:nrow(coor.2d)){

		# INDICES OF VIEWS THAT ARE NOT NA
		v_idx <- which(!is.na(coor.2d[i, seq(1, ncol(coor.2d), 2)])) # index of cameras with non-NA

		if(length(v_idx) < min.views){
			coor.3d[i, ] <- c(NA,  NA,  NA)
			rmse[i] <- NA
			next
		}

		# SET UP MATRICES FOR LINEAR ALGEBRA SOLUTION		
		A <- matrix(0, length(v_idx)*2, 3)
		A[seq(1, length(v_idx)*2, 2), 1] <- coor.2d[i, v_idx*2-1]*cal.coeff[9, v_idx] - cal.coeff[1, v_idx]
		A[seq(1, length(v_idx)*2, 2), 2] <- coor.2d[i, v_idx*2-1]*cal.coeff[10, v_idx] - cal.coeff[2, v_idx]
		A[seq(1, length(v_idx)*2, 2), 3] <- coor.2d[i, v_idx*2-1]*cal.coeff[11, v_idx] - cal.coeff[3, v_idx]
		A[seq(2, length(v_idx)*2, 2), 1] <- coor.2d[i, v_idx*2]*cal.coeff[9, v_idx] - cal.coeff[5, v_idx]
		A[seq(2, length(v_idx)*2, 2), 2] <- coor.2d[i, v_idx*2]*cal.coeff[10, v_idx] - cal.coeff[6, v_idx]
		A[seq(2, length(v_idx)*2, 2), 3] <- coor.2d[i, v_idx*2]*cal.coeff[11, v_idx] - cal.coeff[7, v_idx]

		B <- matrix(0, length(v_idx)*2, 1)
		B[seq(1, length(v_idx)*2, 2), 1] <- cal.coeff[4, v_idx]-coor.2d[i, v_idx*2-1]
		B[seq(2, length(v_idx)*2, 2), 1] <- cal.coeff[8, v_idx]-coor.2d[i, v_idx*2]

		# SOLVE -- INCLUDING EXTRA TRANSFORMATIONS FOR NON-SQUARE MATRIX
		coor.3d[i, ] <- solve(t(A) %*% A) %*% (t(A) %*% B)

		# BACK-CALCULATE 2D COORDINATES FROM RECONSTRUCTED 3D COORDINATES TO ESTIMATE ERROR
		coor_inv <- A %*% as.matrix(coor.3d[i, 1:3])

		# ESTIMATE ROOT MEAN SQUARE RECONSTRUCTION ERROR, WITH PROPER DEGREES OF FREEDOM
		rmse[i] <- sqrt(sum((B - coor_inv)^2) / (length(B)-3))
	}

	l <- list(coor.3d=coor.3d, rmse=rmse)
	class(l) <- 'dltReconstruct'
	l
}

summary.dltReconstruct <- function(object, ...){
	r <- ''
	r <- c(r, '\ndltReconstruct Summary\n')

	r <- c(r, '\tRMSE by Landmark:\n')
	for(i in 1:length(object$rmse)) r <- c(r, '\t\t', rownames(object$coor.3d)[i], ': ', format(object$rmse[i]), '\n')

	class(r) <- "summary.unifyLandmarks"
	r
}

print.summary.dltReconstruct <- function(x, ...) cat(x, sep='')