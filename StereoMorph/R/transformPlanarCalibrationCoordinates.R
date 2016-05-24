transformPlanarCalibrationCoordinates <- function(tpar, nx, ny, sx, sy=NULL){

	# IF Y-SCALE IS NULL, MAKE SAME AS X-SCALE
	if(is.null(sy)) sy <- sx

	# CONVERT INPUT PARAMETERS TO MATRIX
	tpar <- matrix(tpar, ncol=6, byrow=TRUE)

	# MAKE SINGLE SCALED GRID POINT MATRIX, IN REAL-WORLD UNITS
	grid <- cbind(rep(0:(nx-1), ny)*sx, c(matrix(t(matrix((ny-1):0, nrow=ny, ncol=nx)), nrow=1, byrow=F))*sy, rep(0, nx*ny))

	# CENTER GRID ABOUT GRID CENTROID
	grid <- grid - matrix(colMeans(grid), nrow=nrow(grid), ncol=3, byrow=TRUE)

	# 3D TRANSFORMED MATRIX
	coor.3d <- matrix(NA, nrow=nrow(tpar)*nrow(grid), ncol=3)

	for(i in 1:nrow(tpar)){		
		
		# ROTATE GRID
		grid_transform <- grid %*% rotationMatrixZYX_SM(t=c(tpar[i, 1:3]))

		# TRANSLATE GRID
		coor.3d[((i-1)*nrow(grid)+1):(i*nrow(grid)), ] <- grid_transform + matrix(tpar[i, 4:6], nrow=nrow(grid_transform), ncol=3, byrow=TRUE)
	}

	coor.3d
}