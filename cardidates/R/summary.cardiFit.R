`summary.cardiFit` <-
function(object,
  xmin=0, xmax=365, quantile = 0.05, symmetric=FALSE, ...) {
  
  if (!inherits(object, "cardiFit")) stop("use this only with results of fitweibull")
  p   <- object$p
  fit <- object$fit

  ## identify cardinal dates from fitted curves
  smd  <- CDW(object, xmin=xmin, xmax=xmax, quantile=quantile, symmetric=symmetric)

  cat("\n=========== Summary of Cardinal Dates Algorithm 'Weibull' ===========",
   "\n  number of fitted parameters for Weibull function:", length(p),
   "\n  parameter values:", p,
   "\n               r2 =", object$r2,
   "\n  quantile        =", quantile,
   "\n  cardinal dates  (tMid, tBegin, tEnd):", smd$x,
   "\n  original ymax   =", object$ymax,
   "\n\n"
  )
}

