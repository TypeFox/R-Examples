dltCoefficients <- function(coor.3d, coor.2d){
	# Translated and modified from Matlab function dlt_computeCoefficients written by Ty Hedrick
	# Potential additional methods available at: http://kwon3d.com/theory/dlt/dlt.html

	# CONVERT TO 3D ARRAY WITH ONE VIEW IF MATRIX
	if(length(dim(coor.2d)) == 2) coor.2d <- array(coor.2d, dim=c(nrow(coor.2d), ncol(coor.2d), 1))

	# CHECK THAT DIMENSIONS ARE CORRECT
	if(dim(coor.2d)[1] != dim(coor.3d)[1]) stop(paste0("The number of elements in the first dimension of coor.2d (", dim(coor.2d)[1], ") does not match the number of elements in the first dimension of coor.3d (", dim(coor.3d)[1], ")"))

	# DLT CALIBRATION COEFFICIENT MATRIX, 11 PER VIEW
	cal.coeff <- matrix(NA, 11, dim(coor.2d)[3])

	# INTERNAL ROOT MEAN SQUARE ERROR
	rmse <- rep(NA, dim(coor.2d)[3])

	# LOOP THROUGH AND FIND COEFFICIENTS FOR EACH VIEW
	for(i in 1:dim(coor.2d)[3]){
	
		# SKIP IF ALL NA
		if(sum(!is.na(coor.2d[, , i])) == 0) next

		if(is.null(rownames(coor.3d))){

			# IF ROWNAMES ARE ABSENT ASSUME THAT ROWS CORRESPOND, FIND ROWS WHERE NEITHER NA
			non_na <- (is.na(coor.2d[, 1, i]) == FALSE) + (is.na(coor.3d[, 1]) == FALSE) == 2

			# GET ROWS WHERE NEITHER IS NA
			coor_2d <- coor.2d[non_na, , i]
			coor_3d <- coor.3d[non_na, ]
		}else{

			# FIND COMMON POINTS BETWEEN 3D AND 2D COORDINATE SETS
			rownames_common <- rownames(coor.3d)[rownames(coor.3d) %in% rownames(na.omit(coor.2d[, , i]))]
	
			# GET ONLY POINTS PRESENT IN BOTH MATRICES
			coor_2d <- coor.2d[rownames_common, , i]
			coor_3d <- coor.3d[rownames_common, ]
		}
		
		# COEFFICIENT MATRIX
		A <- matrix(0, 2*nrow(coor_3d), 11)
		for(j in 1:nrow(coor_3d)){
			A[2*j-1, 1:3] <- coor_3d[j, 1:3]
			A[2*j, 5:7] <- coor_3d[j, 1:3]
			A[2*j-1, 4] <- 1
			A[2*j, 8] <- 1
			A[2*j-1, 9:11] <- coor_3d[j, 1:3] * -coor_2d[j, 1]
			A[2*j, 9:11] <- coor_3d[j, 1:3] * -coor_2d[j, 2]
		}

		# SOLUTION MATRIX
		B <- c(as.matrix(t(coor_2d)))

		# SOLVE -- INCLUDING EXTRA TRANSFORMATIONS FOR NON-SQUARE MATRIX
		cal.coeff[, i] <- solve(t(A) %*% A) %*% (t(A) %*% B)
		
		# USE CALIBRATION COEFFICIENTS AND 3D COORDINATES TO BACK-CALCULATE PIXEL COORDINATES
		coor_2d_inv <- dltInverse(cal.coeff=cal.coeff[, i], coor.3d=coor_3d)

		# FIND ROOT MEAN SQUARED ERROR BETWEEN COEFFICIENT-DERIVED AND ACTUAL PIXEL COORDINATES
		rmse[i] <- sqrt(mean((coor_2d_inv-coor_2d)^2))
	}

	l <- list(cal.coeff=cal.coeff, rmse=rmse)
	class(l) <- 'dltCoefficients'
	l
}

summary.dltCoefficients <- function(object, ...){

	r <- ''
	r <- c(r, '\ndltCoefficients Summary\n')

	r <- c(r, '\tCalibration RMSE:\n')
	for(i in 1:length(object$rmse)){
		r <- c(r, '\t\t', i, ': ', format(object$rmse[i]), '\n')
	}

	class(r) <- "summary.dltCoefficients"
	r
}

print.summary.dltCoefficients <- function(x, ...) cat(x, sep='')