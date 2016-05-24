plot.synchrony <- function (x, main="", xlab="Values from randomizations", ylab="Frequency", 
                            line.col="red", lty=2, lwd=1, col="grey", ...) {
  if (!is.null(x$rands) & class(x)=="synchrony") {
    if (!is.null(x$w.corrected))
      x$obs=x$w.corrected
    hist(x$rands, main=main, xlab=xlab, ylab=ylab, col=col, ...)
    abline(v=x$obs, col=line.col, lty=lty, lwd=lwd, ...)
    box()
  }
  else {
    stop("No permutation data to plot")
  }
}