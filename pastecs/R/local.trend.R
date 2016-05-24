"local.trend" <-
function(x, k=mean(x), plotit=TRUE, type="l", cols=1:2, ltys=2:1, xlab="Time", ylab="cusum", ...) {
	call <- match.call()
	Data <- deparse(substitute(x))
	if (!is.null(ncol(x)))
		stop("only univariate series are allowed")
	if (length(x) < 3)
		stop("you need at least 3 values in the series")
	x <- as.ts(x)
	x2 <- cumsum(x-k)
	# put x at the same scale as x2
	xmin <- min(x)
	xmax <-max(x)
	x2min <- min(x2)
	x2max <-max(x2)
	x <- (x - xmin) / (xmax - xmin) * (x2max - x2min) + x2min
	x2.ts <- ts(x2, frequency=frequency(x), start=start(x))
	if (plotit == TRUE) {
		if (length(cols) < 2) cols <- rep(cols, 2)
		if (length(ltys) < 2) ltys <- rep(ltys, 2)
		plot(x, type=type, col=cols[1], lty=ltys[1], xlab=xlab, ylab=ylab, ...)
		lines(x2.ts, col=cols[2], lty=ltys[2])
	}
	res <- x2.ts
	attr(res, "k") <- k
	attr(res, "data") <- Data
	attr(res, "call") <- call
	class(res) <- c("local.trend", class(x2.ts))		# turn it into a 'local.trend' object
	res
}
