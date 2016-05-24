plot.variofit <-
  function (x, xlab="Lag distance", ylab="Variogram", 
            col.pts="black", col.line="red", pch=21, ...) {
    
    plot(x$bins, x$vario, xlab=xlab, ylab=ylab, type="p", col=col.pts, pch=pch, ...)
    lines(x$bins, x$fit, col=col.line)
    legend(x="topleft", legend=paste(x$model, " AIC: ", format(x$AIC, dig=2), sep=""), 
           bty="n", lty=1, col=col.line)
  }
