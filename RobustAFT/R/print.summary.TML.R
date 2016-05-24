print.summary.TML<-function(x, digits = max(3, getOption("digits") - 3), 
  signif.stars = getOption("show.signif.stars"), ...)
{
	cat("\nCall:", deparse(x$call), "\n\n")
	if(length(cf <- coef(x))){
		if(nsingular <- sum(x$aliased))
			cat("\nCoefficients:(", nsingular, " not defined because of singularities)\n", sep="")
		else cat("\nCoefficients:\n")
		printCoefmat(cf, digits = digits, signif.stars = signif.stars, na.print = "NA",...)
	}
	else cat("No coefficients\n\n")
	cat("\n")
	cat("Total number of observations: ", length(x$residuals), "\n")
	cat("Number of retained observations: ", x$tn, "\n\n")
	if(length(x$sigma)) cat("Scale Estimate:", x$sigma, "\n")
	else cat("No scale estimate\n\n")
	cat("\n")
	cat("Final cut-off values: ", x$cutoff.values, "\n")
	if(x$call[1]=="TML.censored()")  cat("\nFitted with a robust accelerated failure time model with", sQuote(x$errors), "errors\n")
	if(x$call[1]=="TML.noncensored()")    cat("\nFitted by truncated maximum likelihood with", sQuote(x$errors), "errors\n")
	invisible(x)
}