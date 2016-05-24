"disto" <-
function(x, max.dist=nrow(x)/4, plotit=TRUE, disto.data=NULL) {
	if (is.null(disto.data)) {	# Calculate distogram
		call <- match.call()
		data <- deparse(substitute(x))
		x <- as.matrix(x)
		if (is.null(ncol(x)) || ncol(x) < 2)
			stop("There must be at least two columns (series) in the dataset")
		n <- nrow(x)
		if (is.null(n) || n < 10)
			stop("There must be at least 10 observations in the series")
		max.dist <- round(max.dist)
		if (max.dist < 0) max.dist <- round(n/3)
		if( max.dist >= n) max.dist <- n-1
		distance <- dist(1:n)
		x2 <- x^2
		val <- dist((x2 / apply(x2, 1, sum))^.5)^2
		val <- data.frame(distance=as.numeric(distance), distogram=as.numeric(val))
    	# Calculate mean values for each distance
    	res <- rep(0, max.dist)
    	for (i in 1:max.dist) {
    		res[i] <- mean(val[val$distance == i,]$distogram, na.rm=TRUE)/2	
    	}
    	res <- list(distance=1:max.dist, distogram=res)
    	res <- as.data.frame(res)
    	attr(res, "data") <- data
    	attr(res, "call") <- call
    } else {		# Use disto.data instead
    	res <- disto.data
    }
    if (plotit == TRUE) {	# plot the distogram
    	plot(res$distance, res$distogram, type="l", xlab="distance", ylab="delta^2", main=paste("Distogram for:", attr(res, "data")))
    }
    res
}
