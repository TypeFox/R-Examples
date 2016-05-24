checkCornerOrder <- function(corners, nx, ny, print.progress, verbose){

	if(sum(rowSums(corners) == 0) > 0 || nrow(corners) != nx*ny){
		if(print.progress && verbose) cat(paste0("\torderCorners() failed.\n"))
		return(matrix(0, nrow=1, ncol=2))
	}

	plot <- FALSE
	
	# CHECK THAT ALL POINTS ARE NEAR LINE OF CORRESPONDING ROW
	for(i in seq(1, nx*ny, by=nx)){
		dpl <- distancePointToLine(corners[i:(i+nx-1), ], corners[i, ], corners[i+nx-1, ])
		if(max(dpl) > 20){
			if(print.progress && verbose) cat(paste0("\torderCorners() failed: at least one point more than 20 pixels from expected checkerboard line.\n"))
			return(return(matrix(0, nrow=1, ncol=2)))
		}
	}

	# FLIP CORNERS TO MATCH IMAGE FOR CLARITY IN DEBUGGING
	corners <- cbind(corners[, 2], -corners[, 1])
	
	# MAKE SURE FIRST ROW ADVANCES TOWARD THE RIGHT (POSITIVE X)
	v1 <- corners[nx, ] - corners[1, ]

	# FLIP NX DIRECTION
	if(v1[1] < 0){
		m <- matrix(1:(nx*ny), nrow=ny, ncol=nx, byrow=TRUE)
		new_order <- c(t(m[, ncol(m):1]))
		corners <- corners[new_order, ]
	}
	
	if(plot){
		plot(corners, asp=1, col=rainbow(nrow(corners), end=0.8))
		points(corners, type='l')
	}

	# FLIP CORNERS BACK TO MATCH IMAGE
	cbind(-corners[, 2], corners[, 1])
}