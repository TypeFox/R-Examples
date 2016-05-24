##Script generated in:
# 2011
# 6:29:12 PM
#by: 
# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

centroidClosure <- function (points3D, w = 2) {
	w <- w + .05
	normV <- function (v)	return(sqrt(sum(v ^ 2)))
	findRadius <- function(points3D, G) {
		D <- apply(points3D, 1, function(x) normV(x - G))
		return(max(D))
	}
	
	n <- nrow(points3D)
	G <- colMeans(points3D)
	radius <- findRadius(points3D, G)
	
	#cannot extend less than the trivial (null) extension
	if( w < max( normV(points3D[1,	] - G), normV(points3D[n, ] - G) ) / radius)
		return('the supplied parameter w is not correct. Any k > 1 will do.')
	
	first.point <- G + w * radius * ((points3D[1, ] - G) / 
				(normV(points3D[1,	] - G)))
	last.point <- G + w * radius * ((points3D[n, ] - G) / 
				(normV(points3D[n, ] - G)))
	points3D <- rbind(as.numeric(first.point), 
			points3D, as.numeric(last.point))
	return(points3D)
}
