print.summary.threshpt <-
function(x,...)
{
	thresh=as.character(x$threshold)
	cat("Call:\n")
	print(x$call)
	cat("\nFormula:\n")
	print(x$formula)
	cat("\nCoefficients:\n")
	printCoefmat(x$coefficients, p.value=TRUE, has.Pvalue=TRUE)
	cat("\n")
	cat("Optimum threshold: ")
	cat(thresh) 
	cat("\nDeviance of the model with optimum threshold: ")
	cat(round(x$deviance,3))
	cat("\n")
}
