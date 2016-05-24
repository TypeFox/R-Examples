'biplot.pcoa' <- 
    function(x, Y=NULL, plot.axes = c(1,2), dir.axis1=1, dir.axis2=1, rn=NULL, ...)
# x = output object from function pcoa.R
# Y = optional sites-by-variables data table
# plot.axes = the two axes to be plotted
# rn = an optional vector, length n, of object labels
# dir.axis.1 = -1 to revert axis 1 for the projection of points and variables
# dir.axis.2 = -1 to revert axis 2 for the projection of points and variables
#
# Author: Pierre Legendre, January 2009
	{
	k <- ncol(x$vectors)
	if(k < 2) stop("There is a single eigenvalue. No plot can be produced.")
	if(k < plot.axes[1]) stop("Axis",plot.axes[1],"does not exist.")
	if(k < plot.axes[2]) stop("Axis",plot.axes[2],"does not exist.")

	if(!is.null(rn)) rownames(x$vectors) <- rn
	labels = colnames(x$vectors[,plot.axes])
	diag.dir <- diag(c(dir.axis1,dir.axis2))
	x$vectors[,plot.axes] <- x$vectors[,plot.axes] %*% diag.dir

	if(is.null(Y)) {
		limits <- apply(x$vectors[,plot.axes], 2, range) 
		ran.x <- limits[2,1] - limits[1,1]
		ran.y <- limits[2,2] - limits[1,2]
		xlim <- c((limits[1,1]-ran.x/10), (limits[2,1]+ran.x/5)) 
		ylim <- c((limits[1,2]-ran.y/10), (limits[2,2]+ran.y/10))

		par(mai = c(1.0, 1.0, 1.0, 0.5))
		plot(x$vectors[,plot.axes], xlab=labels[1], ylab=labels[2], xlim=xlim, ylim=ylim, asp=1)
		text(x$vectors[,plot.axes], labels=rownames(x$vectors), pos=4, cex=1, offset=0.5)
		title(main = "PCoA ordination", line=2.5)
	
		} else {
		# Find positions of variables in biplot:
		# construct U from covariance matrix between Y and standardized point vectors
		# (equivalent to PCA scaling 1, since PCoA preserves distances among objects)
		n <- nrow(Y)
		points.stand <- scale(x$vectors[,plot.axes])
		S <- cov(Y, points.stand)
		U <- S %*% diag((x$values$Eigenvalues[plot.axes]/(n-1))^(-0.5))
		colnames(U) <- colnames(x$vectors[,plot.axes])

		par(mai = c(1, 0.5, 1.4, 0))
		biplot(x$vectors[,plot.axes], U, xlab=labels[1], ylab=labels[2])
		title(main = c("PCoA biplot","Response variables projected","as in PCA with scaling 1"), line=4)
	}
    invisible()
}
