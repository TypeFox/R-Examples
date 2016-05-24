# $Id: $

###_ + TDR, matrix method (for .depthFilter)

".plotZOCfilters" <- function(x, zoc.filter, xlim, ylim,
                              ylab="Depth (m)", ...)
{
    ## Value: Nothing; a plot
    ## --------------------------------------------------------------------
    ## Arguments: x=TDR object; zoc.filter=matrix of filters (as returned
    ## by .depthFilter); xlim=numeric vector of length 2 with axis limit
    ## (defaults to time range of input); ylim=numeric vector of length 2
    ## (upper, lower) with axis limits (defaults to range of input);
    ## ...=passed to legend()
    ## --------------------------------------------------------------------
    ## Purpose: Plot to help finding parameters for .depthFilter
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    if (!is(x, "TDR")) stop ("x is not a TDR object")
    time <- getTime(x)
    depth <- getDepth(x)
    if (length(time) != nrow(zoc.filter)) {
        stop ("x and zoc.filter must have the same number of records")
    }
    nfilters <- ncol(zoc.filter)
    if (missing(xlim)) xlim <- range(time)
    xat <- pretty(xlim)
    xlabels <- format(xat)
    if (missing(ylim)) {
        ylim <- rev(-range(depth))
    } else ylim <- rev(-ylim)
    if (nfilters < 2) stop("zoc.filter must have at least 2 columns")
    if (nfilters < 4) {
        npanels <- 3
        filterCols <- 2
    } else {
        npanels <- 3
        filterCols <- c(2, nfilters - 1)
    }
    filter.colors <- seq(3, 4)
    old.par <- par(no.readonly=TRUE)
    on.exit(par(old.par))
    par(mar=c(3, 4, 0, 1) + 0.1, cex=1.1, las=1)
    layout(seq(npanels))
    plot(time, -depth, type="l", col="gray", ylim=ylim, ylab=ylab,
         axes=FALSE)
    axis(1, at=xat, labels=xlabels); axis(2); box()
    abline(h=0, lty=2)
    legend("bottomleft", legend="input", lty=1, col="gray", ...)
    plot(time, -zoc.filter[, 1], type="l", col=2, ylim=ylim, ylab=ylab,
         axes=FALSE)
    axis(1, at=xat, labels=xlabels); axis(2); box()
    abline(h=0, lty=2)
    if (nfilters > 2) {
        for (f in filterCols) {
            lcolor <- filter.colors[filterCols %in% f]
            lines(time, -zoc.filter[, f], col=lcolor)
        }
        legend("bottomleft", legend=colnames(zoc.filter)[c(1, filterCols)],
               lty=1, col=c(2, filter.colors), ...)
    } else {
        legend("bottomleft", legend=colnames(zoc.filter)[1], lty=1,
               col=c(2, filter.colors), ...)
    }
    plot(time, -zoc.filter[, nfilters], type="l", ylim=ylim, ylab=ylab,
         axes=FALSE)
    axis(1, at=xat, labels=xlabels); axis(2); box()
    abline(h=0, lty=2)
    legend("bottomleft",
           legend=paste("input -", colnames(zoc.filter)[nfilters - 1]),
           lty=1, ...)
}

###_ + TDR, TDRcalibrate method (for general comparisons)

".plotZOCtdrs" <- function(x, y, xlim, ylim, ylab="Depth (m)", ...)
{
    ## Value: Nothing; a plot
    ## --------------------------------------------------------------------
    ## Arguments: x=TDR object; y=TDRcalibrate object; xlim=numeric vector
    ## of length 2 with axis limit (defaults to time range of input);
    ## ylim=numeric vector of length 2 (upper, lower) with axis limits
    ## (defaults to range of input); ...=passed to legend()
    ## --------------------------------------------------------------------
    ## Purpose: Plot to help assessing ZOC
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    if (!is(x, "TDR")) stop ("x is not a TDR object")
    if (!is(y, "TDRcalibrate")) stop ("y is not a TDRcalibrate object")
    time.in <- getTime(x)
    depth.in <- getDepth(x)
    time.out <- getTime(getTDR(y))
    depth.out <- getDepth(getTDR(y))
    if (length(time.in) != length(time.out)) {
        stop("x and y must have the same number of records")
    }
    if (missing(xlim)) xlim <- range(time.in)
    xat <- pretty(xlim)
    xlabels <- format(xat)
    if (missing(ylim)) {
        ylim <- rev(-range(depth.in))
    } else ylim <- rev(-ylim)
    plot(time.in, -depth.in, type="l", col="gray", xlim=xlim, ylim=ylim,
         ylab=ylab, axes=FALSE, las=1)
    axis(1, at=xat, labels=xlabels); axis(2); box()
    abline(h=0, lty=2)
    lines(time.out, -depth.out, col=rgb(200, 0, 0, 100, maxColorValue=255))
    abline(h=0, lty=2)
    legend("bottomleft", legend=c("input", "ZOC depth"),
           lty=1, col=c("gray", "red"), ...)
}


###_ + Emacs local variables
## Local variables:
## allout-layout: (+ : 0)
## End:
