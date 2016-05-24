print.coef.plsRmodel <- function(x,...)
{
  if(!is.null(x$Coeffs)){
    cat("Coefficients of the components\n")
    print(x$CoeffC)
    cat("Coefficients of the predictors (original scale)\n")
    print(x$Coeffs)
  }
  if(!is.null(x$Std.Coeffs)){
    cat("Coefficients of the components\n")
    print(x$CoeffC)
    cat("Coefficients of the predictors (scaled scale)\n")
    print(x$Std.Coeffs)
  }
}
