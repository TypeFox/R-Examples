"print.drm" <-
  function (x, ...) 
{
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat(paste("\nResidual deviance:", round(x$deviance, 3), "\tAIC:", 
            round(x$aic, 3), "\n"))
  if (x$code > 2) 
    cat(paste("WARNING: Convergence code", x$code, "; See ?nlm for details\n"))
  if (any(x$fitted.profiles < 0)) 
    cat("WARNING: Negative fitted profile probabilities; model is wrongly specified\n")
  if (any(unlist(x$fitted.conditionals) < 0) || any(unlist(x$fitted.conditionals) > 
                                                    1)) 
    cat("WARNING: Negative fitted conditional probabilities; model is wrongly specified\n")
}

