"print.floated" <-
  function(x, digits=max(3, getOption("digits") - 3), level = 0.95, ...)
{
  K <- qnorm((1+level)/2)
    
  n <- length(x$coef)
  mat <- matrix("", n, 4)
  ci.mat <- matrix(0, n, 2)

  cm <- x$coefmat
  cat("Floating treatment contrasts for factor ", x$factor, "\n\n")
  mat[,1] <- names(x$coef)
  se <- sqrt(x$var)
  ci.mat[, 1] <- x$coef - K * se
  ci.mat[, 2] <- x$coef + K * se
    
  mat[,2] <- format(x$coef, digits=digits)
  mat[,3]   <- format(se, digits=digits)
  ci.mat <- format(ci.mat, digits=digits)
  mat[,4] <- paste("(", ci.mat[,1], ",", ci.mat[,2], ")", sep="")

  dimnames(mat) <- list(rep("", n), c("Level", "Coefficient",
                                      "Std. Error", "95% Floating CI"))
  print(mat, quote=FALSE)
  cat("\nError limits over all contrasts: ", 
      paste(format(c(0.99, x$limits),
                   digits=2)[-1], collapse=","),"\n")
}
