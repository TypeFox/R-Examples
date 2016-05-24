distanceGridUnits <- function(pairs, nx){

	# EMPTY DISTANCE MATRIX
	d <- rep(NA, length=nrow(pairs))

	for(i in 1:nrow(pairs)){
		
		# FIND POSITIONS OF GRID POINT PAIRS
		xy1 <- c((pairs[i, 1] - 1) %% nx, ceiling(pairs[i, 1]/nx) - 1)
		xy2 <- c((pairs[i, 2] - 1) %% nx, ceiling(pairs[i, 2]/nx) - 1)

		# FIND INTERPOINT DISTANCES
		d[i] <- sqrt(sum((xy1 - xy2)^2))
	}

	d
}