.graypr <- function(x.axis=TRUE, y.axis=TRUE, x.major=TRUE, y.major=TRUE, x.minor=TRUE, y.minor=TRUE, x.malty=1, y.malty=1, x.milty=1, y.milty=1){
    if (x.axis)
        axis(1, lwd=0, lwd.ticks=1)
    if (y.axis)
        axis(2, lwd=0, lwd.ticks=1)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], border=NA, col=gray(0.85))
    x.ticks <- axTicks(1)
    y.ticks <- axTicks(2)
    if (x.major){
        abline(v=x.ticks, col=gray(0.90), lty=x.malty, lwd=2)
    }
    if (y.major){
        abline(h=y.ticks, col=gray(0.90), lty=y.malty, lwd=2)
    }
    if (x.minor){
        x.sep <- diff(x.ticks)[1]/2
        x.minorgrid <- c(min(x.ticks)-x.sep, x.ticks+x.sep)
        abline(v=x.minorgrid, col=gray(0.90), lty=x.milty)
    }
    if (y.minor){
        y.sep <- diff(y.ticks)[1]/2
        y.minorgrid <- c(min(y.ticks)-y.sep, y.ticks+y.sep)
        abline(h=y.minorgrid, col=gray(0.90), lty=y.milty)
    }
}
