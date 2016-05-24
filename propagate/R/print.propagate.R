print.propagate <- function(x, ...)
{
  object <- x
  
  ## print error propagation results
  cat("Results from error propagation:\n")
  print(object$prop)
  
  ## print simulation results
  if (!is.na(x$resSIM)) {
    cat("\nResults from Monte Carlo simulation:\n")
    print(object$sim)
  }
}
