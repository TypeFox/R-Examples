print.TransModel <-
function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(x$coefficients)
	cat("\nCovariance Matrix:\n")
	print(x$vcov)
	class(x)<-"TransModel"
}
