"plotAutoC" <-
function(node, plot = TRUE, colour = c("red", "blue", "green", "yellow", "black"), 
    lwd = 5, main = NULL, ...)
#   Plot auto correlation function for single component of OpenBUGS name
{
    sM <- samplesMonitors(node)
    if(length(sM) > 1 || sM != node)
        stop("node must be a scalar variable from the model, for arrays use samplesAutoC")
    nodeName <- sQuote(node)
    sample <- samplesSample(node)
    chain <- samplesGetFirstChain()
    if (sd(sample) > 1.0E-10)
        acfresult <- acf(sample, col = colour[chain], main = if(is.null(main)) nodeName else main, 
            lwd = lwd, demean = TRUE, plot = plot, ...)
    else stop("ACF cannot be computed/plotted: standard deviation <= 1.0E-10")
    acfresult$series <- node
    if(plot) invisible(acfresult)
    else return(acfresult)
}
