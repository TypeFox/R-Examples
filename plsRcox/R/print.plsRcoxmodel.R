print.plsRcoxmodel <- function(x,...)
{
  cat("Number of required components:\n")
  print(x$nt)
  cat("Number of successfully computed components:\n")
  print(x$computed_nt)
  cat("Coefficients:\n")
  print(x$Coeffs)
  cat("Information criteria and Fit statistics:\n")
  print(x$InfCrit)
  if (!is.null(x$family))
  {
    cat("Model with all the required components:\n")
    print(x$FinalModel)
  }
}
