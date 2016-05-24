dltTestCalibration <- function(cal.coeff, coor.2d, nx, sq.size, reciprocal = TRUE){

	# IF SINGLE ASPECT, ADD EXTRA DIMENSION TO MAKE 4D ARRAY
	if(length(dim(coor.2d)) == 3) coor.2d <- array(coor.2d, dim=c(dim(coor.2d)[1:2], 1, dim(coor.2d)[3]))

	# FIND SECOND GRID DIMENSION
	ny <- dim(coor.2d)[1]/nx

	# GET SQUARE SIZES AND UNITS
	sq.size.num <- as.numeric(gsub('[[:alpha:], ]', '', sq.size))
	sq.size.units <- gsub('[[:digit:]., ]', '', sq.size)

	# EMPTY STRUCTURES FOR RESULTS
	ipd_error <- matrix(NA, nrow=floor((nx*ny)/2), ncol=dim(coor.2d)[3], dimnames=list(NULL, dimnames(coor.2d)[[3]]))
	ipd <- matrix(NA, nrow=floor((nx*ny)/2), ncol=dim(coor.2d)[3], dimnames=list(NULL, dimnames(coor.2d)[[3]]))
	epipolar_error <- matrix(NA, nrow=nx*ny, ncol=dim(coor.2d)[3], dimnames=list(NULL, dimnames(coor.2d)[[3]]))
	aitr_dist_error <- epipolar_error
	adj_pair_ipd_error <- matrix(NA, nrow=ny*length(seq(1, nx-1, by=2)), ncol=dim(coor.2d)[3], dimnames=list(NULL, dimnames(coor.2d)[[3]]))
	adj_pair_mean_pos <- array(NA, dim=c(ny*length(seq(1, nx-1, by=2)), 3, dim(coor.2d)[3]), dimnames=list(NULL, c('x','y','z'), dimnames(coor.2d)[[3]]))
	aitr_error <- array(NA, dim=c(nx*ny, 3, dim(coor.2d)[3]), dimnames=list(dimnames(coor.2d)[[1]], c('x','y','z'), dimnames(coor.2d)[[3]]))
	aitr_pos <- aitr_error
	aitr_centroid_dist <- matrix(NA, nrow=nx*ny, ncol=dim(coor.2d)[3], dimnames=list(NULL, dimnames(coor.2d)[[3]]))
	adj_pair_centroid_dist <- matrix(NA, nrow=ny*length(seq(1, nx-1, by=2)), ncol=dim(coor.2d)[3], dimnames=list(NULL, dimnames(coor.2d)[[3]]))
	aitr_rmse <- matrix(NA, nrow=3, ncol=dim(coor.2d)[3], dimnames=list(c('x','y','z'), dimnames(coor.2d)[[3]]))

	# LOOP THROUGH EACH ASPECT
	for(aspect in 1:dim(coor.2d)[3]){

		# GET 3D RECONSTRUCTED POINTS
		dlt_rec <- dltReconstruct(cal.coeff, coor.2d[, , aspect,])


		## FIND ALIGNED TO IDEAL RECONSTRUCTION ERRORS
		# MAKE THEORETICAL GRID OF SAME SIZE FOR ESTIMATE COMPARISON
		coor_3d <- transformPlanarCalibrationCoordinates(tpar=rep(0, 6), nx=nx, ny=ny, sx=sq.size.num)

		# GET OPTIMAL POINT ALIGNMENT
		coor_3d_unify <- findOptimalPointAlignment(dlt_rec$coor.3d, coor_3d)
		
		# SAVE 3D COORDINATE POSITIONS
		aitr_pos[, , aspect] <- coor_3d_unify

		# SAVE ERROR IN REFERENCE-ESTIMATE POINT POSITION AND POSITION OF ESTIMATE POINTS
		aitr_error[, , aspect] <- dlt_rec$coor.3d - coor_3d_unify


		## FIND INTERPOINT DISTANCE ERROR
		# GENERATE RANDOM POINT PAIRS, NO POINTS ARE REPEATED
		ipd_list <- findInterpointDistanceError(coor.3d=dlt_rec$coor.3d, nx=nx, ny=ny, sq.size=sq.size.num)
		ipd[, aspect] <- ipd_list$ipd
		ipd_error[, aspect] <- ipd_list$ipd.error
		adj_pair_ipd_error[, aspect] <- ipd_list$adj.pair.ipd.error
		adj_pair_mean_pos[, , aspect] <- ipd_list$adj.pair.mean.pos


		## FIND EPIPOLAR ERROR
		# MAKE MATRIX FOR PAIRING BETWEEN FIRST AND SUBSEQUENT VIEWS
		ee_mat <- matrix(NA, nrow=nx*ny, ncol=dim(coor.2d)[4]-1)

		# FIND EPIPOLAR ERROR BETWEEN FIRST VIEW AND SUBSEQUENT OTHER VIEWS
		for(k in 2:dim(coor.2d)[4])
			ee_mat[, k-1] <- dltEpipolarDistance(p1=coor.2d[, , aspect, 1], p2=coor.2d[, , aspect, k], cal.coeff[, c(1,k)], reciprocal=reciprocal)

		# ADD TO MATRIX
		epipolar_error[, aspect] <- rowMeans(ee_mat)
	}

	# GET CENTROID OF AITR POINTS
	aitr_centroid <- colMeans(apply(aitr_pos, 2, matrix, byrow=TRUE))

	# GET CENTROID OF AITR POINTS
	adj_pair_mean_centroid <- colMeans(apply(adj_pair_mean_pos, 2, matrix, byrow=TRUE))

	for(aspect in 1:dim(coor.2d)[3]){
	
		# GET DISTANCE FROM CENTROID FOR ESTIMATE POINTS
		aitr_centroid_dist[, aspect] <- distancePointToPoint(aitr_centroid, aitr_pos[, , aspect])

		# GET DISTANCE FROM CENTROID FOR ADJOINING PAIR MEAN POSITIONS
		adj_pair_centroid_dist[, aspect] <- distancePointToPoint(adj_pair_mean_centroid, adj_pair_mean_pos[, , aspect])
		
		# CALCULATE DISTANCE BETWEEN EACH ESTIMATED AND REFERENCE POINT
		aitr_dist_error[, aspect] <- sqrt(rowSums(aitr_error[, , aspect]^2))

		# CALCULATE AITR RMS ERROR FOR EACH ASPECT
		aitr_rmse[, aspect] <- sqrt(colMeans(aitr_error[, , aspect]^2))
	}

	# CALCULATE AITR RMS ERROR FROM ALL TEST ORIENTATIONS
	aitr_dist_rmse <- sqrt(colMeans(aitr_dist_error^2))

	# CALCULATE EPIPOLAR RMS ERROR FROM ALL TEST ORIENTATIONS
	epipolar_rmse <- sqrt(colMeans(epipolar_error)^2)

	# CALCULATE RMS ERROR OF INTERPOINT DISTANCE BASED ON ALL TEST POINT DISTANCES
	ipd_rmse <- sqrt(colMeans(ipd_error^2))

	l <- list(
		num.aspects=dim(coor.2d)[3],
		num.views=dim(coor.2d)[4],
		sq.size.num=sq.size.num,
		sq.size.units=sq.size.units,
		epipolar.error=epipolar_error, 
		epipolar.rmse=epipolar_rmse, 
		ipd.error=ipd_error, 
		pair.dist=ipd,
		ipd.rmse=ipd_rmse, 
		adj.pair.ipd.error=adj_pair_ipd_error, 
		adj.pair.mean.pos=adj_pair_mean_pos, 
		adj.pair.centroid.dist=adj_pair_centroid_dist,
		aitr.error=aitr_error,
		aitr.dist.error=aitr_dist_error,
		aitr.dist.rmse=aitr_dist_rmse,
		aitr.rmse=aitr_rmse,
		aitr.pos=aitr_pos,
		aitr.centroid.dist=aitr_centroid_dist
		)
	class(l) <- 'dltTestCalibration'
	l
}

summary.dltTestCalibration <- function(object, print.tab = '', ...){

	r <- ''

	r <- c(r, '\n', print.tab, 'dltTestCalibration Summary\n')
	
	r <- c(r, print.tab, '\tNumber of aspects: ', object$num.aspects, '\n')
	r <- c(r, print.tab, '\tNumber of views: ', object$num.views, '\n')
	r <- c(r, print.tab, '\tSquare size: ', object$sq.size.num, ' ', object$sq.size.units, '\n')
	r <- c(r, print.tab, '\tNumber of points per aspect: ', nrow(object$epipolar.error), '\n')
	r <- c(r, print.tab, '\tAligned ideal to reconstructed (AITR) point position errors:\n')
	r <- c(r, print.tab, '\t\tAITR RMS Errors (X, Y, Z): ')
	r <- c(r, paste0(paste(format(rowMeans(object$aitr.rmse)), collapse=paste0(' ', object$sq.size.units, ', ')), ' ', object$sq.size.units))
	r <- c(r, '\n')

	r <- c(r, print.tab, '\t\tMean AITR Distance Error: ', format(mean(object$aitr.dist.error)), ' ', object$sq.size.units, '\n')
	r <- c(r, print.tab, '\t\tAITR Distance RMS Error: ', format(mean(object$aitr.dist.rmse)), ' ', object$sq.size.units, '\n')

	r <- c(r, print.tab, '\tInter-point distance (IPD) errors:\n')
	r <- c(r, print.tab, '\t\tIPD RMS Error: ', format(mean(object$ipd.rmse)), ' ', object$sq.size.units, '\n')
	r <- c(r, print.tab, '\t\tIPD Mean Absolute Error: ', format(mean(abs(object$ipd.error))), ' ', object$sq.size.units, '\n')
	r <- c(r, print.tab, '\t\tMean IPD error: ', format(mean(object$ipd.error)), ' ', object$sq.size.units, '\n')

	r <- c(r, print.tab, '\tAdjacent-pair distance errors:\n')
	r <- c(r, print.tab, '\t\tMean adjacent-pair distance error: ', format(mean(object$adj.pair.ipd.error)), ' ', object$sq.size.units, '\n')
	r <- c(r, print.tab, '\t\tMean adjacent-pair absolute distance error: ', format(mean(abs(object$adj.pair.ipd.error))), ' ', object$sq.size.units, '\n')
	r <- c(r, print.tab, '\t\tSD of adjacent-pair distance error: ', format(sd(object$adj.pair.ipd.error)), ' ', object$sq.size.units, '\n')

	r <- c(r, print.tab, '\tEpipolar errors:\n')
	r <- c(r, print.tab, '\t\tEpipolar RMS Error: ', format(mean(object$epipolar.rmse, na.rm=TRUE)), ' px\n')
	r <- c(r, print.tab, '\t\tEpipolar Mean Error: ', format(mean(object$epipolar.error, na.rm=TRUE)), ' px\n')
	r <- c(r, print.tab, '\t\tEpipolar Max Error: ', format(max(object$epipolar.error, na.rm=TRUE)), ' px\n')
	r <- c(r, print.tab, '\t\tSD of Epipolar Error: ', format(sd(object$epipolar.error, na.rm=TRUE)), ' px\n')

	class(r) <- "summary.dltTestCalibration"
	r
}

print.summary.dltTestCalibration <- function(x, ...) cat(x, sep='')