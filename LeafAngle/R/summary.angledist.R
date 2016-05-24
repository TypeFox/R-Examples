`summary.angledist` <-
function(object,...){

	cat("----- Leaf angle distribution -----\n\n")
	cat("Distribution:",object$distribution,"\n")
	cat("Fit with:", object$fitmethod, "\n") 
	if(object$fitmethod == "Chi-squared")cat("Chi-squared statistic =",round(object$chisq,3),"\n")
	if(object$fitmethod == "Log-likelihood")cat("Log-likelihood =",round(object$loglik,3),"\n")
	if(object$fitmethod == "Log-likelihood")cat("AIC =", round(object$AIC,3),"\n")
	if(object$distribution %in% c("ellipsoid","rotatedell"))
			cat("Parameter X =", object$distpars, "\n")
	if(object$distribution == "twoparbeta")
			cat("Parameter tmean =", object$distpars[1], ", parameter tvar =", object$distpars[2], "\n")
}

