# Script comments and history
# 2011
# 5:35:25 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

edgeIntersections <-
		function (points3D, indices, e) 
{
#	if (missing(points3D)) 
#		stop("edgeIntersections: Argument 'points3D' missing, with no default\n")
#	if (missing(indices)) 
#		stop("edgeIntersections: Argument 'indices' missing, with no default\n")
#	if (missing(e)) 
#		stop("edgeIntersections: Argument 'e' missing, with no default\n")
	
	deuc <- function(x,y) 
		sqrt(sum((x - y)^2))
	
	e.max <- length(indices) - 1
	if (e > e.max) 
		stop("edgeIntersections: Argument 'e' out of bound\n")
	e.indices <- indices[e + c(0, 1)]
	set <- indices[-length(indices)]
	int.out <- c()
	for (i in set) {
		pointsij <- points3D[c(e.indices, i, i + 1), ]
		d.1 <- deuc(pointsij[1, 1:2], pointsij[2, 1:2])
		d.2 <- deuc(pointsij[3, 1:2], pointsij[4, 1:2])
		if(deuc(pointsij[1, 1:2], pointsij[3, 1:2]) > (d.1 + d.2))
			int <- 0
		else {
			int <- singleIntersection(pointsij, "sign")
			int.out <- rbind(int.out, c(i, int))
			#colnames(int.out) <- c("edge", "sign")
		}
	}
	return(int.out)
}



edgeIntersectionsK<-
		function (points3D, indices, e) 
{
	#if (missing(points3D)) 
	#	stop("edgeIntersections: Argument 'points3D' missing, with no default\n")
	#if (missing(indices)) 
	#	stop("edgeIntersections: Argument 'indices' missing, with no default\n")
	#if (missing(e)) 
	#	stop("edgeIntersections: Argument 'e' missing, with no default\n")
	
	deuc <- function(x,y) 
		sqrt(sum((x - y)^2))
	
	e.max <- length(indices) - 1
	if (e > e.max) 
		stop("edgeIntersections: Argument 'e' out of bound\n")
	e.indices <- indices[e + c(0, 1)]
	set <- indices[-length(indices)]
	int.out <- c()
	for (i in set) {
		pointsij <- points3D[c(e.indices, i, i + 1), ]
		d.1 <- deuc(pointsij[1, 1:2], pointsij[2, 1:2])
		d.2 <- deuc(pointsij[3, 1:2], pointsij[4, 1:2])
		if(deuc(pointsij[1, 1:2], pointsij[3, 1:2]) > (d.1 + d.2))
			int <- 0
		else {
			int <- singleIntersection(pointsij, "k"); 
			if (length(int)>1) {
				int.out <- rbind(int.out, c(i, int))
				colnames(int.out) <- c("edge", "sign","k");
			}
		}
	}
	return(int.out)
}