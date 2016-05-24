"plot.stat.slide" <-
function(x, stat="mean", col=c(1,2), lty=c(par("lty"), par("lty")), leg=FALSE, llab=c("series", stat), lpos=c(1.5, 10), xlab="time", ylab="y", main=paste("Sliding statistics"), ...) {
	# The next function actually draws the graph
	stat.slide.graph <- function(X, Stat, Col, Lty, Leg, Llab, Lpos, Xlab, Ylab, Main, ...) {
		# Verify if Stat is among possible values
		STATS <- c("min", "max", "median", "mean", "pos.median", "pos.mean", "geo.mean", "pen.mean")
		stat.idx <- pmatch(Stat, STATS)
		if (is.na(stat.idx)) 
			stop("invalid stat value")
		if (stat.idx == -1) 
			stop("ambiguous stat value")
		ysld <- switch(stat.idx,
				"min"=unlist(X$stat["min",]),
				"max"=unlist(X$stat["max",]),
				"median"=unlist(X$stat["median",]),
				"mean"=unlist(X$stat["mean",]),
				"pos.median"=unlist(X$stat["pos.median",]),
				"pos.mean"=unlist(X$stat["pos.mean",]),
				"geo.mean"=unlist(X$stat["geo.mean",]),
				"pen.mean"=unlist(X$stat["pen.mean",]))
		# Verify that there is something in ysld
		if ( sum(is.na(ysld)) == length(ysld))
			stop(paste(Stat, "was not calculated in x!"))
		plot(X$x, X$y, type="l", col=Col[1], lty=Lty[1], xlab=Xlab, ylab=Ylab, main=Main, ...)
		# Construct x and y vectors for the sliding statistics
		xsld <- sort(rep(X$xcut,2))
		yn <- length(ysld)
		ysld[2:(2*yn+1)] <- ysld[floor(seq(1,yn+0.5, by=0.5))]
		ysld[1] <- min(X$x,na.rm=TRUE)
		ysld[2*yn+2] <- min(X$x,na.rm=TRUE)
		lines(xsld, ysld, type="l", col=Col[2], lty=Lty[2], ...)
		# If Leg is TRUE, print a legend
		if (Leg == TRUE) {
			legend(Lpos[1], Lpos[2], Llab, col=Col, lty=Lty)
		}
		
	}
	invisible(stat.slide.graph(x, stat, col, lty, leg, llab, lpos, xlab, ylab, main, ...))
}
