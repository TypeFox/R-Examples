print.smacofJK <- function(x,...)
{
  cat("\nCall: ")
  print(x$call)
  cat("\n")
  cat("SMACOF Jackknife\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("Value loss function:", x$loss, "\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
  cat("Stability measure:",x$stab,"\n")
  cat("Cross validity:",x$cross,"\n")
  cat("Dispersion:",x$disp,"\n")
  cat("\n")
}

