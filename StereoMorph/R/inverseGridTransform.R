inverseGridTransform <- function(grids){

	# CREATE TRANSFORMATION PARAMETER MATRIX - ASSUME GRID IN FIRST POSITION AT 0,0,0
	t.param <- matrix(0, nrow=dim(grids)[3], 6)

	# FIND TRANSLATION PARAMETERS
	centroids <- matrix(NA, nrow=dim(grids)[3], 3)
	
	# FIND CENTROID OF EACH GRID
	for(i in 1:nrow(centroids)) centroids[i, ] <- colMeans(grids[, , i])

	# ADD TRANSLATIONS FROM FIRST TO OTHER CENTROIDS TO TRANSFORMATION PARAMETER MATRIX
	t.param[, 4:6] <- centroids - matrix(centroids[1, ], nrow=nrow(centroids), ncol=ncol(centroids), byrow=TRUE)

	# CENTER EACH GRID AT THE ORIGIN
	grids_c <- grids
	for(i in 1:dim(grids_c)[3]) grids_c[, , i] <- grids[, , i] - matrix(centroids[i, ], nrow=dim(grids)[1], ncol=dim(grids)[2], byrow=TRUE)

	# FOR NOW, JUST RETURN ANGLES ABOUT THE Z,Y AND X AXES
	uvr <- grids_c[1, , 1]/sqrt(sum(grids_c[1, , 1]^2))
	xyr <- uvr * c(1,1,0)
	xzr <- uvr * c(1,0,1)
	yzr <- uvr * c(0,1,1)

	for(i in 2:dim(grids_c)[3]){

		uv <- grids_c[1, , i]/sqrt(sum(grids_c[1, , i]^2))
		xy <- uv * c(1,1,0)
		xz <- uv * c(1,0,1)
		yz <- uv * c(0,1,1)

		t.param[i, 1:3] <- c(avectors(xy, xyr), avectors(xz, xzr), avectors(yz, yzr))
	}

	t.param
}