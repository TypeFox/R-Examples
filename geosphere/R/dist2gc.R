# based on code by Ed Williams
# Licence: GPL
# http://williams.best.vwh.net/avform.htm#XTE
# Port to R by Robert Hijmans
# October 2009
# version 0.1
# license GPL3

dist2gc <- function(p1, p2, p3, r=6378137) {
	toRad <- pi / 180 
	p1 <- .pointsToMatrix(p1)
	p2 <- .pointsToMatrix(p2)
	p3 <- .pointsToMatrix(p3)
	r <- as.vector(r)
	p <- cbind(p1[,1], p1[,2], p2[,1], p2[,2], p3[,1], p3[,2], r)
	p1 <- p[,1:2]
	p2 <- p[,3:4]
	p3 <- p[,5:6]
	r <- p[,7]
	
	tc <- bearing(p1, p2) * toRad
	tcp <- bearing(p1, p3) * toRad
    dp <- distCosine(p1, p3, r=1)
	xtr <- asin(sin(tcp-tc) * sin(dp)) * r
	as.vector(xtr)
}

