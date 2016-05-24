"CenterD2" <-
function(series, window=nrow(series)/5, plotit=TRUE, add=FALSE, type="l", level=0.05, lhorz=TRUE, lcol=2, llty=2,...) {
	call <- match.call()
	data <- deparse(substitute(series))
	if (is.null(class(series)) || class(series)[1] != "mts")
		stop("series must be a multiple regular time series object")
	Unit <- attr(series, "units")
	UnitTxt <- GetUnitText(series)
	# Test the length of the serie, range and step...
	n <- nrow(series)
	if (window[1] < 1 || window[1] > n)
	   	stop("window must be larger or equal to 1, and smaller or equal to n")
	Lags.vec <- 1:n
	D2.vec <- Lags.vec
	# Calculate CenterD2 for each lag
	x <- as.matrix(series)
	w <- scale(x[1:window,])
	R <- solve(cor(w))
	for (i in 1:window) {
		v <- w[i,]
		D2.vec[i] <- (t(v) %*% R) %*% v
	}
	for (j in 1:(n-window)) {
		w <- scale(x[j+(1:window),])
		R <- solve(cor(w))
		v <- w[window,]
		D2.vec[window+j] <- (t(v) %*% R) %*% v
	}
	res <- list(lag=Lags.vec, D2=D2.vec)
	as.data.frame(res)
	res$call <- call
	res$data <- data
	res$type <- "CenterD2"
	res$window <- window
	res$level <- level
	res$chisq <- qchisq(1-level, ncol(series))
	res$units.text <- UnitTxt
	attr(res, "units") <- Unit
		
	# Do we need to plot the graph?
	if (plotit == TRUE) {
		if (add == TRUE) {
			lines(res$lag, res$D2, ...)
		} else {
			plot(res$lag, res$D2, type=type, xlab=paste("lag (", UnitTxt, ")", sep=""), ylab="D2", main=paste("CenterD2 for:", data), ...)
			if (lhorz == TRUE) {
				abline(h=res$chisq, col=lcol, lty=llty)
			}
		}
	}
	class(res) <- "D2"
	res 	# Return results
}
