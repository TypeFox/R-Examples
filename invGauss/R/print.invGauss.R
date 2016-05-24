"print.invGauss" <-
function(x, ...)
{
	cat("Call:\n")
	print(x$call)	#
# coefficients <- cbind(Value = x$coefficients, "Std. Error" = x$SE, Z = x$coefficients/x$SE)
	cat("\nCoefficients:\n")
	print(x$coefficients, ...)
	invisible(x)
}

