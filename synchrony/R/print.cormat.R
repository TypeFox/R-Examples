print.cormat <- function (x, digits=max(3L, getOption("digits") - 3L), ...) {
  cat("Specified mean pairwise correlation: ")
  cat(format(x$rho, digits=digits))
  
  cat("\nSpecified standard deviation: ")
  cat(format(x$sigma, digits=digits))
  
  cat("\nSpecified mean: ")
  cat(format(x$mu, digits=digits))
  
  cat("\nCommunity matrix: \n")
  print(x$community, digits=digits)
}
