distPointToPoint <- function(p1, p2 = NULL){

	# CONVERT VECTORS TO MATRICES
	if(is.vector(p1) && is.null(p2)) p1 <- matrix(p1, nrow=length(p1), ncol=1)
	if(is.vector(p1) && !is.null(p2)) p1 <- matrix(p1, nrow=1, ncol=length(p1))
	if(is.vector(p2)) p2 <- matrix(p2, nrow=1, ncol=length(p2))

	if(is.null(p2)){

		# SINGLE MATRIX OF POINTS WITH MORE THAN TWO ROWS AND ONE COLUMN - FIND INTERPOINT DISTANCES
		if(is.matrix(p1) && nrow(p1) >= 2 && ncol(p1) == 1) return(abs(p1[2:length(p1), ] - p1[1:(length(p1)-1), ]))

		# SINGLE MATRIX OF POINTS WITH MORE THAN TWO ROWS - FIND INTERPOINT DISTANCES
		if(is.matrix(p1) && nrow(p1) >= 2){
			d <- rep(NA, nrow(p1)-1)
			for(i in 2:nrow(p1)) d[i-1] <- distPointToPoint(p1[i-1, ], p1[i, ])
			return(d)
		}

		# SINGLE POINT, RETURN DISTANCE FROM ORIGIN
		return(sqrt(sum(p1^2)))
	}

	# FOR ONE POINT TO MANY, MATCH MATRIX DIMENSIONS
	if(nrow(p1) == 1 & nrow(p2) > 1) p1 <- matrix(p1, nrow=nrow(p2), ncol=ncol(p2), byrow=TRUE)
	if(nrow(p2) == 1 & nrow(p1) > 1) p2 <- matrix(p2, nrow=nrow(p1), ncol=ncol(p1), byrow=TRUE)

	# PAIRWISE DISTANCES BETWEEN POINTS IN TWO MATRICES
	sqrt(rowSums((p1 - p2)^2))
}