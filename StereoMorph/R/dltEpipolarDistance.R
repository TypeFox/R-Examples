dltEpipolarDistance <- function(p1, p2, cal.coeff, reciprocal=FALSE){

	if(reciprocal){		
		# GET EPIPOLAR DISTANCE
		d1 <- dltEpipolarDistance(p1=p1, p2=p2, cal.coeff=cal.coeff)

		# GET RECIPROCAL EPIPOLAR DISTANCE
		d2 <- dltEpipolarDistance(p1=p2, p2=p1, cal.coeff=cal.coeff[, ncol(cal.coeff):1])
		
		# AVERAGE TWO DISTANCE SETS
		return(rowMeans(cbind(d1, d2)))
	}

	# CONVERT INPUTS TO MATRIX IF NOT ALREADY
	if(!is.matrix(p1)){p1 <- matrix(p1, ncol=2)}
	if(!is.matrix(p2)){p2 <- matrix(p2, ncol=2)}
	
	# VECTOR FOR SAVING EPIPOLAR DISTANCES
	d <- rep(NA, length=nrow(p2))

	# IF P1 IS SINGLE POINT, FIND DISTANCE FROM P1 TO ONE OR MANY P2
	if(nrow(p1) == 1){

		# GET EPIPOLAR LINE OF P1 IN VIEW 2
		r <- dltEpipolarLine(p1, cal.coeff)

		# GET DISTANCE FROM POINT TO LINE
		for(i in 1:nrow(p2)) d[i] <- distancePointToLine(p2[i, ], c(0, r$b), c(1, r$m + r$b))
	}

	# IF MULTIPLE P1 AND P2 FIND EPIPOLAR DISTANCE BETWEEN PAIRS
	if(nrow(p1) > 1 && nrow(p1) == nrow(p2)){

		for(i in 1:nrow(p1)){

			# FIND EPIPOLAR LINE OF P1 IN VIEW 2
			r <- dltEpipolarLine(p1[i, ], cal.coeff)

			# GET DISTANCE FROM POINT TO LINE
			d[i] <- distancePointToLine(p2[i, ], c(0, r$b), c(1, r$m + r$b))
		}
	}

	d
}