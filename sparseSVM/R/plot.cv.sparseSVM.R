# plot.cv.sparseSVM

plot.cv.sparseSVM <- function(x, log.l = TRUE, nvars = TRUE, ...) {
  l <- x$lambda
  if (log.l) {
    l <- log(l)
    xlab <- expression(log(lambda))
  } else {
    xlab <- expression(lambda)
  }
  
  ## Calculate y
  L.cve <- x$cve - x$cvse
  U.cve <- x$cve + x$cvse
  y <- x$cve
  L <- L.cve
  U <- U.cve
  
  if (x$eval.metric == 'me') {
    ylab <- "Cross-validation Prediction Error"
  } else {
    ylab <- "Cross-validation Error"
  }

  ylim <- range(c(L, U))
  ind <- ((U-L)/diff(ylim) > 1e-3)
  plot.args = list(x=l, y=y, ylim=ylim, xlab=xlab, ylab=ylab, type="n", xlim=rev(range(l)), las=1)
  new.args = list(...)
  if (length(new.args)) {
    plot.args[names(new.args)] <- new.args
  }
  do.call("plot", plot.args)
  suppressWarnings(arrows(x0=l[ind], x1=l[ind], y0=L[ind], y1=U[ind], 
                          code=3, angle=90, col="gray80", length=.03))
  points(l, y, col="red", pch=19, cex=.5)
  if (nvars) {
    n.s <- apply(coef(x$fit, lambda=x$lambda)!=0, 2, sum)-1
    axis(3, at=l, labels=n.s, tick=FALSE, line=-0.5)
    mtext("Variables selected", cex=0.8, line=1.5)
  }
}
