"plotDensity" <-
function(node, plot=TRUE, main = NULL, xlab = "" , ylab = "", col = "red", ...)
#   Plot posterior density for single component of OpenBUGS name
{
    sM <- samplesMonitors(node)
    if(length(sM) > 1 || sM != node)
        stop("node must be a scalar variable from the model, for arrays use samplesDensity")
    nodeName <- sQuote(node)
    sampleSize <- samplesSize(node)
    sample <- samplesSample(node)

    absSample <- abs(sample)
    intSample <- as.integer(absSample + 1.0E-10)
    zero <- absSample - intSample
    intSample <- as.integer(sample)
    if (sum(zero) > 0){
        if (is.R())
          d <- density(sample, adjust = 1.25)
        else
          d <- density(sample)
        if (plot)
            plot(d$x, d$y, type = "l", main = if(is.null(main)) nodeName else main,
                 xlab = xlab , ylab = ylab, col = col, ...)
        res <- d
    }
    else{
        histogram <- table(intSample) / sampleSize
        xRange <- range(intSample)
        xLim <- c(xRange[1] - 0.5, xRange[2] + 0.5)
        if (plot)
            plot(histogram, type = "h", xlim = xLim, ylim = c(0, 1),
                 main = if(is.null(main)) nodeName else main,
                 xlab = xlab , ylab = ylab, col = col, ...)
        res <- histogram
    }
    if (plot) invisible(res) else return(res)
}
