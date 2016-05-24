`plot.cardiFit` <-
function(x, y=NULL,
  xmin=0, xmax=365, quantile = 0.05, symmetric=FALSE, ...) {
  
  if (!inherits(x, "cardiFit")) stop("use this only with results of fitweibull")
  p   <- x$p
  fit <- x$fit
  f   <- x$ymax

  ## select the adequate version of fweibull
  if (length(p) == 4) {
    fweibull <- fweibull4
  } else {
    fweibull <- fweibull6
  }

  ## identify cardinal dates from fitted curves
  smd  <- CDW(x, xmin=xmin, xmax=xmax, quantile=quantile, symmetric=symmetric)

  ## plot data, curve and cardinal dates
  plot(fit$x, fit$y * f, xlab="x values", ylab="backtransformed y values", ...)
  xnew <- seq(min(fit$x), max(fit$x), length=100)
  lines(xnew, fweibull(xnew, p) * f, col="darkgreen", lwd=2)
  abline(v=smd$x, col="grey", lty="dashed")
  points(smd$x, smd$y * f, col="tomato", pch=16, cex=1.4)
}

