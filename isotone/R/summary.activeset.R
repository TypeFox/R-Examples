summary.activeset <- function(object, ...)
{
  cat("\nCall:\n")
  print(object$call)
  cat("\nLoss value:", object$fval, "\n")
  cat("Number of iterations: ", object$niter, "\n")

  cat("\nLangrange Multipliers:\n")
  print(round(object$lambda, 3))
  
  if (!is.null(object$isocheck))
  {
    cat("\nKKT feasibility check:\n")
    chdf <- as.data.frame(unlist(object$isocheck))
    colnames(chdf) <- c("Value")
    print(round(chdf, 3))
  }
  cat("\n")
}