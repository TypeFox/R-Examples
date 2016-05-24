print.varComp <-
function(x, ...)
{
	cat("Variance component model fit", '\n')
	cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
	if(inherits(x$fixef, 'varCompFixEf')) {
		print(x$fixef)
	}else{
		cat("\nFixed effect estimates:", '\n')
		print(coef(x, 'fixed'))
	}
	cat("\nVariance component estimates:", '\n')
	print(coef(x, 'varComp'))
	cat(sprintf("\nNumber of observations: %d\n", nobs(x)))
	invisible(x)
}

summary.varComp = function(object, ...)
{## FIXME: add more summary information
	object$fixef = fixef(object,...)
	class(object) = c('summary.varComp', class(object))
	object
}
