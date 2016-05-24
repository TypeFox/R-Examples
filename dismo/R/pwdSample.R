# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : April 2011
# Version  1.0
# Licence GPL v3


pwdSample <- function(fixed, sample, reference, tr=0.33, nearest=TRUE, n=1, lonlat=TRUE, warn=TRUE) {

	distHaversine <- function (p1, p2) {
		r <- 6378137
		toRad <- pi/180
		p1 <- p1 * toRad
		p2 <- p2 * toRad
		p <- cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2])
		dLat <- (p[, 4] - p[, 2])
		dLon <- (p[, 3] - p[, 1])
		a <- sin(dLat/2) * sin(dLat/2) + cos(p[, 2]) * cos(p[, 4]) * sin(dLon/2) * sin(dLon/2)
		dist <- 2 * atan2(sqrt(a), sqrt(1 - a)) * r
		as.vector(dist)
	}

	distGeo <- function (x, y) {
		n <- nrow(x)
		m <- nrow(y)
		dm <- matrix(ncol = m, nrow = n)
		for (i in 1:n) {
			dm[i, ] <- distHaversine(x[i, ,drop=FALSE], y)
		}
		return(dm)
	}
	
	
	distPlane <- function (x, y) {
		dfun <- function(x, y) {
			sqrt( (x[,1] - y[,1])^2 + (x[,2] - y[,2])^2 )
		}
		n = nrow(x)
		m = nrow(y)
		dm = matrix(ncol = m, nrow = n)
		for (i in 1:n) {
			dm[i, ] = dfun(x[i, ,drop=FALSE], y)
		}
		return(dm)
	}
	
	if (lonlat) {
		distfun <- distGeo
	} else {
		distfun <- distPlane
	}
	
	stopifnot(tr > 0)
	
	n <- round(n)
	stopifnot(n >= 1)
	
	if (inherits(fixed, 'SpatialPoints')) fixed <- coordinates(fixed)
	if (inherits(sample, 'SpatialPoints')) sample <- coordinates(sample)
	if (inherits(reference, 'SpatialPoints')) reference <- coordinates(reference)
	fixed     <- as.matrix(fixed)[,1:2]
	sample    <- as.matrix(sample)[,1:2]
	reference <- as.matrix(reference)[,1:2]

	if (warn) {
		if (nrow(sample) < nrow(fixed)) {
			warning("nrow(sample) < nrow(fixed)")
		}
	}
	
	fromd <- apply(distfun(fixed, reference), 1, min)
	tod <- apply(distfun(sample, reference), 1, min)

	ngb <- matrix(NA, nrow=length(fromd), ncol=n)

	iter <- sample(1:nrow(fixed)) # random order
	
	for (j in 1:n) {
		for (i in iter) {
			d <- abs(tod - fromd[i])
			if (min(d) < (tr * fromd[i])) {
				if (nearest) {
					x <- which.min(d)
				} else {
					x <- sample(which(d < (tr * fromd[i])), size=1) 
				}
				ngb[i, j] <- x
				tod[x] <- Inf
			} 
		}
	}
	return(ngb)
}


