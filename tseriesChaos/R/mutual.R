#Author: Antonio, Fabio Di Narzo. Last Modified $Date: 2005-12-18 16:27:25 +0100 (dom, 18 dic 2005) $
mutual <- function(series, partitions=16, lag.max=20, plot=TRUE, ...) {
	series <- (series-min(series))/(diff(range(series)))
	corr <- numeric(lag.max+1)
	for(i in 0:lag.max) {
		hist <- matrix(0, partitions, partitions)
		hist <- .C("mutual", series=as.double(series), length=as.integer(length(series)), lag=as.integer(i), partitions=as.integer(partitions), hist=as.double(hist), PACKAGE="tseriesChaos")[["hist"]]
		hist <- matrix(hist, partitions, partitions)/sum(hist)
		histx <- apply(hist, 1, sum)
		hist <- hist[hist!=0]
		histx<- histx[histx!=0]
		corr[i+1] <- sum(hist*log(hist)) - 2*sum(histx*log(histx))
	}
	names(corr) <- paste(0:lag.max)
	class(corr) <- "ami"
	if(plot) plot.ami(corr, ...)
	invisible(corr)
}

plot.ami <- function(x, main=NULL, ...) {
	if(is.null(main)) main <- "Average Mutual Information"
	plot(0:(length(x)-1), x, type="h", xlab="lag", ylab="AMI", main=main)
	abline(h=0)
}
