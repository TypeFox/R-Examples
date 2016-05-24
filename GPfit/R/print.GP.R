## A Print command for a class GP opbject
##
## May 8th, 2012

print.GP <- function(x, ...){
if (is.GP(x) == FALSE){
	stop("The object in question is not of class \"GP\" \n")
}
corr = x$correlation_param

cat("\nNumber Of Observations: n = ")
cat(nrow(x$X))
cat("\nInput Dimensions: d = ")
cat(ncol(x$X))
cat("\n\n")
if (corr$type == "exponential"){
	cat("Correlation: Exponential (power = ", corr$power, ")",sep="")
} else {
	cat("Correlation: Matern (nu = ", corr$nu, ")",sep="")
}
cat("\n")
cat("Correlation Parameters: \n")
beta_val = data.frame(beta_hat = t(x$beta), row.names = '[1]')
print(beta_val, ...)
cat("\n")
cat("sigma^2_hat: ")
print(x$sig2, ...)
cat("\n")
cat("delta_lb(beta_hat): ")
new_delta = c(x$delta)
print(new_delta, ...)
cat("\n")
cat("nugget threshold parameter: ")
cat(x$nugget_threshold_parameter)
cat("\n\n")
}