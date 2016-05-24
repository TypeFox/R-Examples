#' Scale vectors to fill the current plot

#' @param vectors a two-column matrix giving the end points of a collection of vectors
#' @param bbox the bounding box of the containing plot region within which the vectors are to be plotted
#' @param origin origin of the vectors
#' @param factor maximum length of the rescaled vectors
#' @return scale factor 

vecscale <- function(vectors, 
		bbox=matrix(par("usr"), 2, 2),
		origin=c(0, 0), factor=0.95) {	
	scale <- c(sapply(bbox[,1] - origin[1], function(dist) dist/vectors[,1]), 
			sapply(bbox[,2] - origin[2], function(dist) dist/vectors[,2])) 
	scale <- factor * min(scale[scale > 0])
	scale
}


