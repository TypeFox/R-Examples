`print.smacofSP` <-
function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("Model:",x$model,"\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("\nStress-1 value:",round(x$stress,3),"\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}

