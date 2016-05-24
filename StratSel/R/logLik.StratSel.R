logLik.StratSel <- function(object, ...) #UseMethod("logLik")
#print.logLik <- function(x, digits = getOption("digits"), ...)
{
	if(missing(digits)) digits = getOption("digits")
	x <- object
	ll <- x$logLik
    cat("'log Lik.' ", paste(format(c(ll), digits = digits), collapse = ", "),
        " (df=", length(coef(x)), ")\n", sep = "")
    invisible(x)
}