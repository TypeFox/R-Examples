#' @rdname plot_imprecise_object
#' @title Plotting Imprecise Objects
#' @description Generic function for plotting \code{imprecise} objects.   
#' @param x an object for which a plot is needed.
#' @param after.obs logical indicating imprecise prior or posterior. 
#' @param ... additional arguments affecting the plot produced. 
#' @export 
plot.impinf <- function(x, after.obs=FALSE, ...){
	
	if ( !after.obs ) vtx <- x$vtx
	else vtx <- x$vtx1

	# 2 dimension
	if ( ncol(vtx) == 2 ){
		graphics::plot(x=vtx, ...)
		graphics::polygon(x=vtx, col="azure2", border="darkblue", lty="dashed", lwd=3, ...)
		graphics::points(x=vtx, col="red", pch=19, cex=0.5)
		graphics::text(x=vtx, labels=rownames(vtx), pos=4)
	}
	
	# 3 dimension requires rgl package 
	if ( ncol(vtx) == 3 ) {
		vtx0 <- x$vtx0
		imat <- x$imat
		suppressWarnings(rgl::triangles3d(x=vtx0[imat,1],y=vtx0[imat,2],z=vtx0[imat,3], color="green"))
		rgl::decorate3d()
		rgl::texts3d(x=vtx[,1],y=vtx[,2],z=vtx[,3], text=rownames(vtx), col="red")
	}
	
}

#' @rdname plot_imprecise_object
#' @title Constructing Probability Band 
#' @param min lower limit of a distribution
#' @param max upper limit of a distribution
#' @export
pbox <- function(x, min, max, ...){
	
	stopifnot(inherits(x, "impinf"))
	vtx1 <- x$vtx1
	ztrunc <- x$ztrunc
	y <- x$y
	
	if(ncol(vtx1)==2) vtx1 <- cbind(0, vtx1)
	graphics::plot(c(0,0), xlim=c(min, max), ylim=c(0,1), type="n", ylab="CDF", xlab=expression(theta), ...) 
	graphics::abline(h=c(0,1), col="grey", lty="dashed")
	t <- seq(from=min, to=max, by=0.2)
	cols <- grDevices::rainbow(nrow(vtx1))
	for(i in 1:nrow(vtx1)){
	 	pars1 <- as.vector(vtx1[i,])
	 	if(!ztrunc) p <- pcpm(q=t, pars=pars1)
	 	if(ztrunc) p <- pcpm.ztrunc(q=t, pars=pars1, ny=length(y))
	 	graphics::lines(t,p, col=cols[i])
	}
}
