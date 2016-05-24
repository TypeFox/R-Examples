"print.depratio" <-
  function (x, digits=max(3, getOption("digits") - 3), ...) 
{
  ord <- x$call$ord
  cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  cat("Observed adjacent ", if (is.null(ord)) 2
  else ord, "-way dependence ratios:\n", sep = "")
  print.default(format(x$tau, digits=digits), print.gap=2,
                quote = FALSE,...)
  if(!is.null(x$call$n.boot) && x$call$n.boot>0){
  cat("\nBootstrap lower confidence levels:\n")
  print.default(format(x$boot$lcl, digits = digits), print.gap = 2, 
                quote = FALSE,...)
   cat("\nBootstrap upper confidence levels:\n")
  print.default(format(x$boot$ucl, digits = digits), print.gap = 2, 
                quote = FALSE,...)
}
}
