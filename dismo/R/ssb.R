# Author: Robert J. Hijmans, r.hijmans@gmail.com
# Date : August 2011
# Version  1.0
# Licence GPL v3


ssb <- function(p, a, reference, lonlat=TRUE, avg=TRUE) {

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
	
	if (inherits(p, 'SpatialPoints')) p <- coordinates(p)
	if (inherits(a, 'SpatialPoints')) a <- coordinates(a)
	if (inherits(reference, 'SpatialPoints')) reference <- coordinates(reference)
	p <- as.matrix(p)[,1:2]
	a <- as.matrix(a)[,1:2]
	reference <- as.matrix(reference)[,1:2]
		
	pdist <- distfun(p, reference)
	adist <- distfun(a, reference)

	pd <- apply(pdist, 1, min)
	ad <- apply(adist, 1, min)
	
	if (avg) {
		res <- cbind(mean(pd), mean(ad))
		colnames(res) <- c('p', 'a')
		return(res)
	} else {
		return( list(pd, ad) )
	}
}




 
 