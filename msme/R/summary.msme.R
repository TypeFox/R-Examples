summary.msme <- function(object, ...) {

### Create a coefficient table  
  z <- with(object, coefficients / se.beta.hat)
  zTable <-
      with(object, 
           data.frame(Estimate = coefficients,
                      SE = se.beta.hat,
                      Z = z,
                      p = 2 * pnorm(-abs(z)),
                      LCL = coefficients - 1.96 * se.beta.hat,
                      UCL = coefficients + 1.96 * se.beta.hat))
  rownames(zTable) <- colnames(object$X)

### Prepare part of the coefficient table for printing  
  z.print <- zTable
  z.print$p <- formatC(z.print$p, digits = 3, format="g")

### Build a list of output objects
  summ <- list(call = object$call,
               coefficients = zTable,
               deviance = object$deviance,
               null.deviance = object$null.deviance,
               df.residual = object$df.residual,
               df.null = object$df.null) 
  if(!is.null(object$pearson)) {
    summ$null.pearson <- object$null.pearson
    summ$pearson <- object$pearson
    summ$dispersion <- object$dispersion
  }
### Write out a set of results
  cat("\nCall:\n")
  print(object$call)
  cat("\nDeviance Residuals:\n")
  print(summary(as.numeric(residuals(object))))
  if(!is.null(object$presiduals)) {
    cat("\nPearson Residuals:\n")
    print(summary(as.numeric(object$presiduals)))    
  }
  cat("\nCoefficients (all in linear predictor):\n")
  print(z.print, digits = 3, ...)
  cat("\nNull deviance:", summ$null.deviance,
      " on ", summ$df.null, "d.f.") 
  cat("\nResidual deviance:", summ$deviance,
      " on ", summ$df.residual, "d.f.") 
  if(!is.null(object$pearson)) {
    cat("\nNull Pearson:", summ$null.pearson,
        " on ", summ$df.null, "d.f.")
    cat("\nResidual Pearson:", summ$pearson,
        " on ", summ$df.residual, "d.f.")
    cat("\nDispersion:", summ$dispersion)
  }
  cat("\nAIC: ", object$aic)
  cat("\n\nNumber of optimizer iterations: ", object$i, "\n\n")

### Return the list but do not print it.  
  return(invisible(summ))
}
