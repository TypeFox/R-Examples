findInterpointDistanceError <- function(coor.3d, nx, ny, sq.size){

	set.seed(42)
	p_r <- 1:(nx*ny)
	p1 <- sample(p_r, floor((nx*ny)/2), replace=F)
	p2 <- sample(p_r[-p1], floor((nx*ny)/2), replace=F)
	pairs <- matrix(NA, nrow=floor((nx*ny)/2), ncol=2)
	pairs <- cbind(p1, p2)
	
	# GET ADJOINING POINTS FOR FIRST ROW
	x1 <- seq(1, nx-1, by=2)
	x2 <- seq(2, nx, by=2)
	
	# EXPAND TO NUMBER OF ROWS AND CUMULATIVE INDEX WITH EACH ROW
	x1_m <- matrix(x1, nrow=ny, ncol=length(x1), byrow=TRUE) + nx*(0:(ny-1))
	x2_m <- matrix(x2, nrow=ny, ncol=length(x2), byrow=TRUE) + nx*(0:(ny-1))

	# GET TWO-COLUMN MATRIX OF PAIRS OF ADJOINING POINTS
	pairs_adjoin <- cbind(matrix(t(x1_m), ncol=1), matrix(t(x2_m), ncol=1))

	# FIND DISTANCE BETWEEN EACH PAIR OF POINTS, SCALED TO KNOWN GRID SIZE IN REAL-WORLD UNITS
	ipd <- distanceGridUnits(pairs, nx)*sq.size
	ipd_adjoin <- distanceGridUnits(pairs_adjoin, nx)*sq.size

	# FIND ERROR IN INTERPOINT DISTANCE
	ipd_pos <- matrix(NA, nrow(pairs), ncol=3)
	ipd_error <- rep(NA, nrow(pairs))
	adj_pair_ipd_error <- rep(NA, nrow(pairs_adjoin))
	adj_pair_mean_pos <- matrix(NA, nrow(pairs_adjoin), ncol=3)

	for(j in 1:nrow(pairs)){
		ipd_error[j] <- distancePointToPoint(coor.3d[pairs[j, 1], ], coor.3d[pairs[j, 2], ]) - ipd[j]
		ipd_pos[j, ] <- colMeans(coor.3d[pairs[j, 1:2], ])
	}

	for(j in 1:nrow(pairs_adjoin))
		adj_pair_ipd_error[j] <- distancePointToPoint(coor.3d[pairs_adjoin[j, 1], ], coor.3d[pairs_adjoin[j, 2], ]) - ipd_adjoin[j]

	# GET MEAN 3D-RECONSTRUCTED POSITION OF ADJOINING PAIRS
	for(j in 1:nrow(pairs_adjoin))
		adj_pair_mean_pos[j, ] <- colMeans(coor.3d[pairs_adjoin[j, 1:2], ])

	list(
		'ipd'=ipd,
		'ipd.pos'=ipd_pos,
		'ipd.error'=ipd_error,
		'adj.pair.ipd.error'=adj_pair_ipd_error,
		'adj.pair.mean.pos'=adj_pair_mean_pos
	)
}