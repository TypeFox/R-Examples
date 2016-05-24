plot.coef.plsRmodel <- function(x,...)
{
  if(!is.null(x$Coeffs)){
    cat("Coefficients of the components\n")
    return(x$CoeffC)
    cat("Coefficients of the predictors (original scale)\n")
    return(x$Coeffs)
  }
  if(!is.null(x$Std.Coeffs)){
    cat("Coefficients of the components\n")
    return(x$CoeffC)
    cat("Coefficients of the predictors (scaled scale)\n")
    return(x$Std.Coeffs)
  }
}
