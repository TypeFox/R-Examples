distPointToPlane <- function(p, n, q){
	# http://mathworld.wolfram.com/Point-PlaneDistance.html
	# q is a point in the plane

	if(is.matrix(p)){
		if(nrow(p) > 1){
			d <- rep(NA, nrow(p))
			for(i in 1:nrow(p)){
				d[i] <- sum((n*-(q - p[i, ])))  / sqrt(sum(n^2))
			}
			return(d)
		}
	}

	sum((n*-(q - p)))  / sqrt(sum(n^2))
}