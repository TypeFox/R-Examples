plot.lacf <-
function (x, plotcor = TRUE, type = "line", lags = 0:min(as.integer(10 * 
    log10(nrow(x$lacf))), ncol(x$lacf) - 1), tcex = 1, lcol = 1, 
    llty = 1, the.time = NULL, plot.it=TRUE, xlab, ylab, ...) 
{

    if (missing(xlab))
	xlab <- "Time"

    nlags <- length(lags)
    ntime <- nrow(x$lacf)
    if (max(lags) + 1 > ncol(x$lacf)) 
        stop("Maximum lag is too high")
    if (length(lcol) == 1) 
        lcol <- rep(lcol, length(lags))
    if (length(llty) == 1) 
        llty <- rep(llty, length(lags))
    if (length(lcol) != length(lags)) 
        stop("Length of lcol vector has to be 1 or the same as the length of the lags vector")
    if (length(llty) != length(lags)) 
        stop("Length of llty vector has to be 1 or the same as the length of the lags vector")
    if (type == "line") {
        if (plotcor == TRUE) {
	    if (plot.it==TRUE)	{

		if (missing(ylab))
			ylab <- "Autocorrelation"
            	plot(c(1, max(ntime)), c(-1, 1), type = "n", xlab = xlab, 
            		ylab = ylab, ...)

            	for (i in 1:nlags) {
                	lines(1:ntime, x$lacr[, 1 + lags[i]], col = lcol[i], 
                  	lty = llty[i])

			pp <- seq(from = 1, to = ntime, length = 5)
			text(pp, x$lacr[pp, 1 + lags[i]], labels = lags[i], 
			  cex = tcex)
			}
		}
	ans <- x$lacr[, 1+ lags]
	dimnames(ans) <- list(NULL, as.character(lags))

	return(invisible(ans))
        }
        else	{
	    # Plot line plot of autocovariance
	    # Get y axis limits
	    #
	    yl <- range(x$lacf[,1+lags])
	    if (plot.it==TRUE)	{

		if (missing(ylab))
			ylab <- "Autocovariance"

            	plot(c(1, max(ntime)), c(yl[1], yl[2]), type = "n", xlab = xlab, 
            		ylab = ylab, ...)

            	for (i in 1:nlags) {
                	lines(1:ntime, x$lacf[, 1 + lags[i]], col = lcol[i], 
                  	lty = llty[i])

			pp <- seq(from = 1, to = ntime, length = 5)
			text(pp, x$lacf[pp, 1 + lags[i]], labels = lags[i], 
			  cex = tcex)
			}
		}
	ans <- x$lacf[, 1+ lags]
	dimnames(ans) <- list(NULL, as.character(lags))

	return(invisible(ans))
	} 
    }
    else if (type == "persp") {
        if (plotcor == TRUE) {
            m <- x$lacr[, lags + 1]
            zlab <- "ACF"
        }
        else {
            m <- x$lacf[, lags + 1]
            zlab <- "ACF (cov)"
        }
	if (plot.it==TRUE)	{

		if (missing(ylab))
			ylab <- "Lag"
		
		persp(x = 1:ntime, y = lags, z = m[, lags + 1], xlab = xlab, 
		    ylab = ylab, zlab = zlab, ...)
		}
	ans <- m
	dimnames(ans) <- list(NULL, lags)
	return(ans)
    }
    else if (type == "acf") {
        if (is.null(the.time)) 
            stop("You have to specify a time point at which you wish to see the autocovariance/correlation. Specify the.time")
        if (plotcor == TRUE) {
	    
            acfvals <- x$lacr[the.time, lags + 1]
	    if (missing(ylab))
		    ylab <- "ACF"
        }
        else {
            acfvals <- x$lacf[the.time, lags + 1]
	    if (missing(ylab))
		    ylab <- "ACF (cov)"
        }
	if (plot.it==TRUE)	{

		if (missing(xlab))
			xlab <- "Lag"
		plot(c(0, max(lags)), c(min(acfvals, 0), 1), type = "n", 
		    xlab = xlab, ylab = ylab, ...)
		segments(x0 = lags, y0 = 0, x1 = lags, y1 = acfvals)
		abline(h = 0)
		}
        return(invisible(acfvals))
    }
}
