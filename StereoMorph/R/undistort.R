undistort <- function(x, image.size, center=c(0,0), k=rep(0, 3), p=rep(0, 2)){

	# FIND SCALE FOR X
	#scale <- 1/max(x)
	
	#x <- scale*x
	if(length(k) < 3) k <- c(k, 0)
	
	k[is.na(k)] <- 0
	p[is.na(p)] <- 0

	xu <- x
	if(length(dim(x)) == 4){
		for(i in 1:dim(x)[4]) xu[, , , i] <- undistort(x=x[, , , i], image.size, center=center, k=k, p=p)
		return(xu)
	}
	if(length(dim(x)) == 3){
		for(i in 1:dim(x)[3]) xu[, , i] <- undistort(x=x[, , i], image.size, center=center, k=k, p=p)
		return(xu)
	}
	
	# Distance between x coordinates and center
	x_center_dist <- sqrt(rowSums((x - matrix(center, nrow=nrow(x), ncol=2, byrow=TRUE))^2))

	# Distance from center to corner
	x_corner_dist <- sqrt(sum((image.size - center)^2))

	# Normalize distances to distance from center to corner
	r <- x_center_dist / x_corner_dist

	for(i in 1:nrow(x)){
		radd <- (1 + k[1]*r[i]^2 + k[2]*r[i]^4 + k[3]*r[i]^6)*(x[i, ]-center)
		tand <- c(
			2*p[1]*(x[i, 2]-center[2]) + p[2]*(r[i]^2 + 2*(x[i, 1]-center[1])^2),
			p[1]*(r[i]^2 + 2*(x[i, 2]-center[2])^2) + 2*p[2]*(x[i, 1]-center[1])
		)
		xu[i, ] <- center + radd + tand
		#xu[i, ] <- center + radd
	}

	#xu <- xu*(1/scale)
	
	xu
}