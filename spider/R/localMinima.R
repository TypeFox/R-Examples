localMinima <- function(distobj){
	den <- density(distobj)
	a <- rep(NA, length(den$y)-2)
	for(i in 2:(length(den$y)-1)) a[i-1] <- den$y[i-1] > den$y[i] & den$y[i+1] > den$y[i]
	den$localMinima <- den$x[which(a)]
	den$data.name <- deparse(substitute(distobj))
	den$call <- paste("density.default(", den$data.name, ")", sep="")
	print(den$localMinima)
	invisible(den)
}