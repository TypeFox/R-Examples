"lines.stat.slide" <-
function(x, stat="mean", col=3, lty=1, ...) {
	# The next function actually draws the lines
	stat.slide.lines <- function(X, Stat, Col, Lty, ...) {
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
		# Construct x and y vectors for the sliding statistics
		xsld <- sort(rep(X$xcut,2))
		yn <- length(ysld)
		ysld[2:(2*yn+1)] <- ysld[floor(seq(1,yn+0.5, by=0.5))]
		ysld[1] <- min(X$x,na.rm=TRUE)
		ysld[2*yn+2] <- min(X$x,na.rm=TRUE)
		lines(xsld, ysld, type="l", col=Col, lty=Lty, ...)
	}
	invisible(stat.slide.lines(x, stat, col, lty, ...))
}
