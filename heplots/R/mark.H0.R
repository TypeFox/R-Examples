# last modified 2 Jan 2010

# mark the location of a point corresponding to a null hypothesis

mark.H0 <- function(x=0, y=0, z=NULL, label, cex=2, pch=19, col="green3", lty=2, pos=2) {
	if (is.null(z)) {
		points(x,y, cex=cex, col=col, pch=pch)
		if (missing(label)) label<-expression(H[0])
		text(x,y, label, col=col, pos=pos)
		if (lty>0) abline(h=y, col=col, lty=lty)
		if (lty>0) abline(v=x, col=col, lty=lty)
	}
	else {
		bbox <- matrix(rgl::par3d("bbox"), nrow=2)
		ranges <- apply(bbox, 2, diff)
		rgl::points3d(x, y, z, size=5*cex, color=col)
		if(lty>0) cross3d(c(x,y,z), (ranges/2), col=col, lty=lty)
	}
}
