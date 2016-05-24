distPointToLine <- function(pt, l1, l2){ # Finds shortest distance between pt and line, not line segment (defined by l1 and l2)

	if(is.matrix(pt)){
		if(nrow(pt) > 1){
			d <- rep(NA, nrow(pt))
			for(i in 1:length(d)){
				d[i] <- sqrt(sum(cprod(pt[i, ] - l1, pt[i, ] - l2)^2)) / sqrt(sum((l2 - l1)^2))
			}
			return(d)
		}
	}

	return(sqrt(sum(cprod(pt - l1, pt - l2)^2)) / sqrt(sum((l2 - l1)^2)))
}