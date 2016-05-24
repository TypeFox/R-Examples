#' Bin plotting
#' 
#' Do a 'scatterplot bin smoothing'
#' 
#' @param x Forcing variable
#' @param y Output
#' @param h the bandwidth (defaults to \code{2*sd(runvar)*length(runvar)^(-.5)})
#' @param nbins number of Bins
#' @param cutpoint Cutpoint
#' @param plot Logical. Whether to plot or only returned silently
#' @param type Whether returns the y averages, or the x frequencies
#' @param xlim,cex,main,xlab,ylab Usual parameters passed to plot(), see \code{\link{par}}
#' @param \ldots further arguments passed to plot. 
#' @return Returns silently values
#' @references McCrary, Justin. 


plotBin <- function(x, y, h = 0.05, nbins = NULL, cutpoint = 0, plot = TRUE, type = c("value", "number"), xlim = range(x, na.rm = TRUE), 
    cex = 0.9, main = NULL, xlab, ylab, ...) {
    
    type <- match.arg(type)
    x_name <- if (missing(xlab)) 
        deparse(substitute(x)) else xlab
    y_name <- if (missing(ylab)) 
        deparse(substitute(y)) else ylab
    
    
    ## Set intervals and midpoints
    min_x <- min(xlim)
    max_x <- max(xlim)
    
    if (!is.null(nbins)) 
        h <- diff(xlim)/nbins
    
    K0 <- ceiling((cutpoint - min_x)/h)  # Number of cells on left
    K1 <- ceiling((cutpoint + max_x)/h)  # Number of cells on right
    K <- K0 + K1
    if (!is.null(nbins) && K != nbins) {
        ranges <- c(cutpoint - min_x, cutpoint + max_x)
        if (which.min(ranges) == 1) {
            K0 <- K0 - 1
        } else {
            K1 <- K1 - 1
        }
        K <- K0 + K1
    }
    
    b_k <- cutpoint - (K0 - c(1:K) + 1) * h  # Lee and Lemieux (2010) p. 308
    mid_points_bk <- b_k + h/2
    n_bins <- length(mid_points_bk)
    brk <- c(b_k, cutpoint + (K1 + 2) * h)
    
    ## compute output (mean of count)
    intervs <- cut(x, breaks = brk, include.lowest = TRUE)
    table_intervs <- table(intervs)
    n_non0_intervs <- sum(table_intervs != 0)
    
    y2 <- switch(type, value = tapply(y, intervs, mean, na.rm = TRUE), number = table_intervs)
    
    
    ## plot
    if (plot) {
        plot(mid_points_bk, as.numeric(y2), pch = 19, cex = cex, xlab = x_name, ylab = y_name, xlim = xlim, ...)
        title(main = main, sub = paste("h=", round(h, 4), ",\\tn bins=", n_non0_intervs, sep = ""))
        abline(v = cutpoint, lty = 2)
    }
    
    ## return invisible result
    res <- data.frame(x = mid_points_bk, y = y2)
    invisible(res)
}



 
