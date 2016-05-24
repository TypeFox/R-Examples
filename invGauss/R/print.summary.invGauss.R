print.summary.invGauss <- function(x, covariance = FALSE, ...){
##
## PRINT A SUMMARIZED invGauss OBJECT
##
cat("Call:\n")
print(x$call)
cat("\nCoefficients:\n")
print(x$coefficients)
cat("\nLog likelihood: ", x$loglik, "\n\n", sep = "")
cat("AIC: ", x$AIC, "\n\n", sep = "")
###	cat("\nIntercept value for centered covariates:\n")
###	print(as.numeric(x$centered.intercept))
###	cat("\nIndividual distances from point of absorption at time 0:\n")
###	print(summary(x$initial.distance))
if(covariance) {
	cat("Asymptotic covariance of coefficients:\n")
	print(x$cov.unscaled)
}
return(invisible(x))
}
