print.TML<-function(x, digits = max(3, getOption("digits") - 3),...)
{
	cat("\nCall:", deparse(x$call), "\n\n")
	if(length(x$th0)){
		cat("Initial Coefficients:\n")
		print.default(format(x$th0, digits = digits), print.gap = 2, quote=FALSE)
	}
	else cat("No initial coefficient\n\n")
	if(length(x$th1)){
		cat("\nFinal Coefficients:\n")
		print.default(format(x$th1, digits = digits), print.gap = 2, quote=FALSE)
	}
	else cat("\nNo final coefficient\n\n")
	if(length(x$v0)) cat("\nInitial Scale: ", x$v0, "\n")
	else cat("\nNo initial scale\n\n")
	if(length(x$v1)) cat("Final Scale: ", x$v1, "\n")
	else cat("No final scale\n\n")
	cat("Final cut-off values: ", x$tl, x$tu, "\n")
	cat("Number of retained observations: ", x$tn, "\n")
	if(x$call[1]=="TML.censored()")  cat("\nFitted with a robust accelerated failure time model with", sQuote(x$errors), "errors\n")
	if(x$call[1]=="TML.noncensored()")    cat("\nFitted by truncated maximum likelihood with", sQuote(x$errors), "errors\n")
	invisible(x)
}