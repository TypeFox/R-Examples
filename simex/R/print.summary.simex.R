print.summary.simex <-
  function (x, digits = max(3, getOption("digits") - 3), ...)
  {
    cat("Call:\n")
    print(x$call)
    cat("\nNaive model: \n")
    print(x$naive.model)
    cat("\nSimex variable :\n")
    
#**Heidi**#    
    if ( length(x$measurement.error) == length(x$SIMEXvariable) ) {
      var.error <- t(matrix(x$measurement.error))
      dimnames(var.error) <- list("Measurement error :", x$SIMEXvariable)
      print(var.error)
    } else {
      m.e <- as.matrix(x$measurement.error)
      var.error <- matrix( ncol = NCOL(x$measurement.error) )
      for ( i in 1:NCOL(x$measurement.error) ) {
        ifelse ( length(unique(m.e[ , i])) == 1, var.error[ , i] <- unique(m.e[ , i]), var.error[ , i] <- "heteroscedastic" )
      }
#      var.error <- t(matrix(rep("heteroscedastic", times = length(x$SIMEXvariable))))
      dimnames(var.error) <- list("Measurement error :", x$SIMEXvariable)
      print(var.error)
    }
#********#    
    
    
    cat("\n\nNumber of iterations: ", x$B, "\n")
    cat("\nResiduals: \n")
    print(summary(x$residuals), digits)
    cat("\nCoefficients: \n")
    if (any(names(x$coefficients) == "asymptotic")) {
      cat("\nAsymptotic variance: \n")
      printCoefmat(x$coefficients$asymptotic, digits = digits)
    }
    if (any(names(x$coefficients) == "jackknife")) {
      cat("\nJackknife variance: \n")
      printCoefmat(x$coefficients$jackknife, digits = digits)
    }
    return(invisible(x))
  }

