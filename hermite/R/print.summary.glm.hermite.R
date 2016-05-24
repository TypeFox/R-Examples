print.summary.glm.hermite <-
  function(x, ...)
  {
    cat("Call:\n")
    print(x$call)
    
    cat("\nDeviance Residuals:\n")
    resid <- c(min(x$resid), quantile(x$resid,c(0.25,0.5,0.75)),max(x$resid))
    names(resid) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(resid)    
    
    cat("\nCoefficients:\n")
    print(x$coefficients)
    cat("(Likelihood ratio test against Poisson is reported by *z value* for *dispersion.index*)\n")
    
    cat("\nAIC: ", x$aic, "\n")
  }
