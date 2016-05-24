mds.plot <- function(X, labels=NULL, col="gray") 
{
	warning("Use of mds.plot(...) is deprecated!")
	return(mdsPlot(X, labels, col));
}

mdsPlot <- function(X, labels=NULL, col="gray") {
	
	if (is.null(labels)) {
		labels <- colnames(X$data)
	}
	
	mds <- cmdscale( X$D, eig=TRUE, k=2  )
	x <- mds$points[,1] 
	y <- mds$points[,2]
	plot(x,y, xlab="Coordinate #1", ylab="Coordinate #2", type="n",
		 main="Multidimensional Scaling")
	
	# determine convex hulls of clusters
	for (cn in unique(labels)) {
		rng <- which(labels==cn)
		idx <- chull(x[rng],y[rng])
		polygon(mds$points[c(rng[idx],rng[idx[1]]),], col=col[which(unique(labels)==cn)])
	
	}
	
	# add labels
	text(x,y, labels)
	
	invisible()
}