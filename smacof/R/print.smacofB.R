`print.smacofB` <-
function(x, ...)
{
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  cat("Model:",x$model,"\n")
  cat("Number of objects:",x$nobj,"\n")
  cat("Stress-1 value:", round(x$stress, 3),"\n")
  cat("Number of iterations:",x$niter,"\n")
  cat("\n")
}

