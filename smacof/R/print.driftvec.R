print.driftvec <- function(x, ...)
{
  cat("\nCall:\n")
  print(x$call, ...)
  cat("\n")
  cat("Symmetric MDS fit:\n")
  cat("Number of objects:",x$fitsym$nobj,"\n")
  cat("Stress-1 value:", round(x$fitsym$stress, 3),"\n")
  cat("Number of iterations:",x$fitsym$niter,"\n")
  cat("\n")
}

