print.summary.fit.linERR <-
function(x, ...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nParameter Summary Table:\n")
    print(x$coefficients)
    
    cat("\nAIC: ", x$aic, "\n")
    cat("Deviance: ", x$dev, "\n")
    cat("Informative risk sets: ", x$inf.rsets, "\n")
  }
