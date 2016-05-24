print.plsRmodel <- function(x,...)
{
  cat("Number of required components:\n")
  print(x$nt)
  cat("Number of successfully computed components:\n")
  print(x$computed_nt)
  cat("Coefficients:\n")
  print(x$Coeffs)
  if (x$typeVC=="none")
  {
    cat("Information criteria and Fit statistics:\n")
    print(x$InfCrit)
  }
  if (x$typeVC %in% c("standard","missingdata","adaptative"))
  {
    cat("Leave one out cross validated PRESS, Information criteria and Fit statistics:\n")
    print(x$InfCrit)
  }
}
