# draw labeled vectors in 3D

vectors3d <- function(x, origin=c(0,0,0), labels=rownames(x), 
	scale=1, 
	col="blue", 
	lwd=1, 
	cex=1,
	 ...) {

	x <- scale*x
#	if (is.vector(origin)) origin <- matrix(origin, ncol=3)

	for(i in 1:nrow(x)) {
		rgl::lines3d( 
			c(origin[1], x[i, 1]),
			c(origin[2], x[i, 2]),
			c(origin[3], x[i, 1]), col=col, lwd=lwd, ...
			)
	}
	if (!is.null(labels)) {
#		if(missing(pos)) pos <- ifelse(x[,1]>0, 4, 2)
  	rgl::texts3d( x, texts=labels, col=col, cex=cex, ...)
		}
}

