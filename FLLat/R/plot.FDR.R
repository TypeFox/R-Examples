plot.FDR <- function(x,xlab="Threshold",ylab="FDR",...) {

    if (!inherits(x,"FDR")) {
        stop("'x' must be of class 'FDR'")
    }

    plot(x=x$thresh.values,y=x$FDRs,xlab=xlab,ylab=ylab,...)

}
