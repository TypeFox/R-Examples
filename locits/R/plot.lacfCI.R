plot.lacfCI <-
function (x, plotcor = TRUE, type = "line", lags = 0:as.integer(10 * 
    log10(nrow(x$lacf))), tcex = 1, lcol = 1, llty = 1, 
    ylim = NULL, segwid = 1, segandcross = TRUE, conf.level = 0.95, 
    plot.it=TRUE, xlab, ylab, sub, ...) 
{
    if (conf.level < 0 || conf.level > 1) 
        stop("conf.level has to be between 0 and 1")
    siz <- 1 - conf.level
    qval <- qnorm(1 - siz/2)
    XCI <- x 
    x <- x$the.lacf
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

		if (missing(xlab))
			xlab <- "Time"

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
        }
	
        else	{
	    # Plot line plot of autocovariance
	    # Get y axis limits
	    #
	    yl <- range(x$lacf[,1+lags])
	    if (plot.it==TRUE)	{


		if (missing(xlab))
			xlab <- "Time"

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
            zlab <- "Autocorrelation"
        }
        else {
            m <- x$lacf[, lags + 1]
            zlab <- "Autocovariance"
        }
	if (plot.it==TRUE)	{
	    if (missing(xlab))
		xlab <- "Time"
	    if (missing(ylab))
		ylab <- "Lag"

            persp(x = 1:ntime, y = lags, z = m[, lags + 1], xlab = xlab, 
                ylab = ylab, zlab = zlab, ...)
	}
    }
    else if (type == "acf") {
	the.time <- XCI$nz
	if (missing(sub))
		sub <- paste("c(", the.time, ", lag)")
	
        if (plotcor == TRUE) {
            acfvals <- x$lacr[the.time, lags + 1]
	    if (missing(ylab))
		    ylab <- "Autocorrelation"
        }
        else {
            acfvals <- x$lacf[the.time, lags + 1]
	    if (missing(ylab))
		    ylab <- "Autocovariance"
        }
        vlags <- XCI$lag
        acvvar <- XCI$cvar
        sv <- match(vlags, lags)
        sw <- 0.2
        x0v <- x1v <- yuv <- ylv <- NULL
        for (i in 1:length(vlags)) {
            if (!is.null(sv[i])) {
                x0v <- c(x0v, vlags[i] - sw/2)
                x1v <- c(x1v, vlags[i] + sw/2)
                yuv <- c(yuv, x$lacf[the.time, vlags[i] + 1] + 
                  qval * sqrt(acvvar[i]))
                ylv <- c(ylv, x$lacf[the.time, vlags[i] + 1] - 
                  qval * sqrt(acvvar[i]))
            }
            else {
                x0v <- c(x0v, NULL)
                x1v <- c(x1v, NULL)
                yuv <- c(yuv, NULL)
                ylv <- c(ylv, NULL)
            }
        }
        if (is.null(ylim)) {
            if (plotcor == FALSE) {
                ylim <- range(c(yuv, ylv, min(acfvals, 0)))
            }
            else ylim <- range(min(acfvals, 0), 1)
        }
        if (plot.it==TRUE)	{
	    if (missing(xlab))
		xlab <- "Lag"
            plot(c(0, max(lags)), c(min(acfvals, 0), 1), type = "n", 
                xlab = xlab, ylab = ylab, ylim = ylim, sub=sub, ...)
            segments(x0 = lags, y0 = 0, x1 = lags, y1 = acfvals, 
                lwd = segwid)
            abline(h = 0)
            if (segandcross == TRUE) 
                points(lags, acfvals, pch = 18)
	
            if (plotcor == FALSE) {
                for (i in 1:length(vlags)) {
                    if (!is.null(sv[i])) {
                      polygon(x = c(x0v[i], x1v[i], x1v[i], x0v[i]), 
                        y = c(ylv[i], ylv[i], yuv[i], yuv[i]), density = 50, 
                        col = rgb(red = 0.9, green = 0.6, blue = 0.6))
                    }
                }
            }
	}
        return(invisible(acfvals))
    }
}
