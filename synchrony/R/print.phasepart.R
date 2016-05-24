print.phasepart <- function (x, digits=max(3L, getOption("digits") - 3L), ...) {
  
  cat("Cross-correlation: ")
  cat(format(x$rho, digits=digits))

  cat("\nAutocorrelation: ")
  cat(format(x$gamma, digits=digits))
  
  cat("\nStandard deviation: ")
  cat(format(x$sigma, digits=digits))
  
  cat("\nMean: ")
  cat(format(x$mu, digits=digits))
  
  #cat("\nTimeseries:\n")
  #print(format(x$timeseries, digits=digits)) 
}
