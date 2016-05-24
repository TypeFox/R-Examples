"print.summary.drm" <-
  function (x, ..., digits = x$digits, correl = FALSE) 
{
  if (is.null(digits)) 
    digits <- options()$digits
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coef)
  cat(paste("\nResidual deviance:", round(x$deviance, digits), 
            "\tAIC:", round(x$aic, digits), "\n\n"))
  cat(paste("Convergence code", x$code, "in", x$niter, "iterations (See help(nlm) for details)\n\n"))
  if (!is.null(x$correl)) {
    correl <- x$correl
    dimnames(correl) <- list(dimnames(x$coefficients)[[1]], 
                             dimnames(x$coefficients)[[1]])
    p <- dim(correl)[2]
    if (p > 1) {
      cat("Correlation of Coefficients:\n")
      correl[!lower.tri(correl)] <- NA
      print(round(correl[-1, -p, drop = FALSE], digits), na = "")
    }
  }
  if (x$code > 2) 
    cat(paste("WARNING: Convergence code", x$code, "; See ?nlm for details\n"))
  if (any(x$fitted.profiles < 0)) 
    cat("WARNING: Negative fitted profile probabilities; model is wrongly specified\n")
  if (any(unlist(x$fitted.conditionals) < 0) || any(unlist(x$fitted.conditionals) > 
                                                    1)) 
    cat("WARNING: Negative fitted conditional probabilities; model is wrongly specified\n")
  invisible()
}

