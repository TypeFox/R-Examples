# Plot methods to add variable vectors showing the 
# original variables in PCA/SVD space.

# Thx: Uwe Ligges for the code for calculating scale...

# for ridge objects, default to first 2 variables
# and show PCA vectors in variable space....
biplot.ridge <-
		function(x, variables=1:2, xlab, ylab, ...) {
	x$svd.V <- t(x$svd.V)
	vnames <- colnames(coef(x))[variables]
	if(missing(xlab)) xlab=vnames[1]
	if(missing(ylab)) ylab=vnames[2]
	
	biplot.pcaridge(x, variables, xlab=xlab, ylab=ylab, ...)
}

biplot.pcaridge <- function(x, variables=(p-1):p, labels=NULL, asp=1, 
		origin, scale, 
		var.lab=rownames(V), 
		var.lwd=1, var.col="black", var.cex=1,
		xlab, ylab,      # override prefix/suffix?
		prefix = "Dim ", # prefix for labels of PCA dimensions
		suffix = TRUE,   # add label suffix with PCA % ?
		...) {
	
	# more convenient versions of arrows() and text()
	Arrows <- function(xy, lenxy, length, angle, col, lwd=1) {
		arrows(xy[1], xy[2], xy[1]+lenxy[,1], xy[2]+lenxy[,2], length=length, angle=angle, lwd=lwd, col=col)
	}
	Text <- function(xy, lenxy, text, col="black", cex=1) {
		text(xy[1]+lenxy[,1], xy[2]+lenxy[,2], text, col=col, cex=cex)
	}
	
	coef <- coef(x)
	p <- ncol(coef)
	if(is.null(x$svd.V)) stop("x must have an svd.V component")
	V <- x$svd.V[,variables]
	
	# add
	pct <- 100*x$svd.D^2 /(sum(x$svd.D^2))
	if (is.logical(suffix) & suffix)
		suffix <- paste( " (", round(pct[variables],3), "%)", sep="" ) else suffix <- NULL
	dimlab <- paste(prefix, variables, suffix, sep="")
	if (missing(xlab)) xlab=dimlab[1]
	if (missing(ylab)) ylab=dimlab[2]
	
	plot(x, variables=variables, labels=labels, asp=asp, xlab=xlab, ylab=ylab, ...)
	
	bbox <- matrix(par("usr"), 2, 2, dimnames=list(c("min", "max"),c("x", "y")))
	if(missing(origin)) origin <- colMeans(bbox)
	
	# plot variable vectors
	if(missing(scale)) {
		scale <- c(sapply(bbox[,"x"] - origin[1], function(dist) dist/V[,1]),
				sapply(bbox[,"y"] - origin[2], function(dist) dist/V[,2]))
		scale <- 0.95* min(scale[scale > 0])
		cat("Vector scale factor set to ", scale, "\n")
	}
	
	Arrows(origin, scale*V, angle=8, length=.1, col=var.col, lwd=var.lwd)
	Text(origin, 1.01*scale*V, var.lab, col=var.col, cex=var.cex)
}

