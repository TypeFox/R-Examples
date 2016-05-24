intersectSphereLine <- function(c, r, x, l, point.compare=NULL){
	# http://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

	# MAKE SURE l IS A UNIT VECTOR
	l <- uvector(l)
	
	# FIND A AND B
	a <- -sum(l*(x - c))
	b <- sum(l*(x - c))^2 - sum((x - c)^2) + r^2
	
	if(b < 0 && abs(b) > 10^-10) stop("No solution possible for sliding joint position.")

	# IF B IS ZERO
	if(abs(b) <= 10^-10) b <- 0

	# DISTANCE ON LINE TO POINT FROM X (LINE ORIGIN)
	d <- c(a + sqrt(b), a - sqrt(b))
	
	# POSSIBLE POINTS
	p <- matrix(NA, nrow=2, ncol=3)
	p[1, ] <- x + d[1]*l
	p[2, ] <- x + d[2]*l

	# FIND DISTANCE FROM TWO INTERSECTION POINTS TO LINE ORIGIN OR COMPARE POINT, IF PROVIDED
	if(is.null(point.compare)){
		dist <- distPointToPoint(x, p)
	}else{
		dist <- distPointToPoint(point.compare, p)
	}

	p[which.min(dist), ]
}