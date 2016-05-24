"vario" <-
function(x, max.dist=length(x)/3, plotit=TRUE, vario.data=NULL) {
	if (is.null(vario.data)) {	# Calculate variogram
		call <- match.call()
		data <- deparse(substitute(x))
		x <- as.ts(x)
		n <- length(x)
		if (n < 10) # Need at least 10 data
			stop("There must be at least 10 observations in the series")
		max.dist <- round(max.dist)
		if (max.dist < 0) max.dist <- round(n/3)
		if( max.dist >= n) max.dist <- n-1
		distance <- dist(1:n)
		val <- outer(x, x, function(X, Y) ((X - Y)^2)/2)
    	val <- val[lower.tri(val)]
    	val <- data.frame(distance=as.numeric(distance), semivario=val)
    	# Calculate mean values for each distance
    	res <- rep(0, max.dist)
    	for (i in 1:max.dist) {
    		res[i] <- mean(val[val$distance == i,]$semivario, na.rm=TRUE)	
    	}
    	res <- list(distance=1:max.dist, semivario=res)
    	res <- as.data.frame(res)
    	attr(res, "data") <- data
    	attr(res, "call") <- call
    } else {		# Use vario.data instead
    	res <- vario.data
    }
    if (plotit == TRUE) {	# plot the variogram
    	plot(res$distance, res$semivario, type="l", xlab="distance", ylab="gamma", main=paste("Semi-variogram for:", attr(res, "data")))
    }
    res
}
