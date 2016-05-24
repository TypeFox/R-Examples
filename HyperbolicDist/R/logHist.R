### Draws a log-histogram
### DJS and Richard Trendall
logHist <- function (x, breaks = "Sturges",
                     include.lowest = TRUE, right = TRUE, 
                     main = paste("Log-Histogram of", xName), 
                     xlim = range(breaks), ylim = NULL, xlab = xName, 
                     ylab = "Log-density", nclass = NULL,
                     htype = "b", ...)
{ 
    xName <- paste(deparse(substitute(x), 500), collapse = "\n")
    histInfor <- hist.default(x, plot=FALSE)
    logDensity <- log(histInfor$density)
    breaks <- histInfor$breaks 
    mids <- histInfor$mids
    nB <- length(breaks)
    counts <- histInfor$counts
    height <- range(logDensity,finite=TRUE)[2] - 
              range(logDensity,finite=TRUE)[1] 
    base <- min(logDensity[is.finite(logDensity)]) - 0.25 *height 

    ## yMax is the max value of logDensity plus another 25%. 
    yMax <- 0.25*abs(max(logDensity))+ max(logDensity) 
    ## plot the log-histogram 
    if(is.null(ylim)) ylim <- range(base,yMax)
    plot(mids, logDensity, xlim = xlim, ylim = ylim,
    type="n", xlab=xlab, ylab=ylab, main=main, ...)
    if (htype=="b"||htype=="p"){
      points(mids, logDensity, ...)
    }
    heights <- rep(0,nB) 

    for (j in  2:(nB-1)) { 
        if(is.finite(max(logDensity[j-1],logDensity[j]))){ 
        heights[j] <- max(logDensity[j-1],logDensity[j]) 
    } 
    else { 
        heights[j] <- NA } 
    } 
    heights[1] <- ifelse(is.finite(logDensity[1]),logDensity[1],NA) 
    heights[nB] <- ifelse(is.finite(logDensity[nB-1]),logDensity[nB-1],NA)
    
    if (htype=="b"||htype=="h"){
      i <- 1:(nB) 
      segments(breaks[i],logDensity[i],breaks[i+1],logDensity[i]) 
      segments(breaks[i],heights[i],breaks[i],base,lty=2) 
      segments(breaks[nB],heights[nB],breaks[nB],base,lty=2)
    }

    r<-list(breaks = breaks, counts = counts, 
        logDensity = logDensity, mids = mids, 
        xName = xName, heights = heights, ylim = ylim) 
    invisible(r) 
} ## End of logHist()
