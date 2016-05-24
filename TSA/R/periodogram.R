`periodogram` <-
function (y, log = "no", plot = TRUE, ylab = "Periodogram", xlab = "Frequency", 
    lwd = 2, ...) 
{
    if (is.matrix(y) && (dim(y)[2] > 1)) 
        stop("y must be a univariate time series")
    sp = spec(y, log = log, plot = FALSE)
    sp$spec = 2 * sp$spec
    temp=sp$spec[sp$freq==.5]
    sp$spec[sp$freq==.5]=temp/2
    if (plot == TRUE) 
        plot(y = sp$spec, x = sp$freq, type = "h", ylab = ylab, 
            xlab = xlab, lwd = lwd, ...)
    return(invisible(sp))
}

