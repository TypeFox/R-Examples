print.scalreg <-
function(x,...)
{
	cat("Call:\n")
	print(x$call)

	if(x$type=="regression"){
		cat("\nEstimated variance:\n")
		print(x$hsigma)
		cat("\nCoefficients:\n")
		print(x$coefficients)
	}
	
	if(x$type=="precision matrix"){
		cat("\nDiagonal:\n")
		print(x$hsigma)
	}
}
