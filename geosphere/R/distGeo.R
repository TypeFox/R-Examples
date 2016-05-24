# Author: Robert J. Hijmans
# Date :  May 2015
# Licence GPL v3


distGeo <- function(p1, p2, a=6378137, f=1/298.257223563) {
	p1 <- .pointsToMatrix(p1) 
	p2 <- .pointsToMatrix(p2) 
	p  <- cbind(p1[,1], p1[,2], p2[,1], p2[,2])
	r <- .Call("inversegeodesic", as.double(p[,1]), as.double(p[,2]), as.double(p[,3]), as.double(p[,4]), as.double(a), as.double(f), PACKAGE='geosphere')
	r <- matrix(r, ncol=3, byrow=TRUE)
	r[, 1]
}

