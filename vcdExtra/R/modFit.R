##
## One-line summary of model fit for a glm/loglm object
##

`modFit` <-
function(x, ...) UseMethod("modFit")

modFit.glm <- function(x, stats="chisq", digits=2, ...) {
	if (!inherits(x,"glm")) stop("modFit requires a glm object")
	result <- NULL
	if ("chisq" %in% stats)
		result <- paste("G^2(",x$df.residual,")=",
				formatC(x$deviance,digits=digits,format="f"),sep="")
	if ("aic" %in% stats)
		result <- paste(result, " AIC=", formatC(x$aic,digits=digits,format="f"),sep="")
	result
}


modFit.loglm <- function(x, stats="chisq", digits=2, ...) {
	if (!inherits(x,"loglm")) stop("modFit requires a loglm object")
	result <- NULL
	if ("chisq" %in% stats)
		result <- paste("G^2(",x$df,")=",
				formatC(x$deviance,digits=digits,format="f"),sep="")
	if ("aic" %in% stats) {
	    aic<-x$deviance-x$df*2
		result <- paste(result, " AIC=", formatC(aic,digits=digits,format="f"),sep="")
		}
	result
}

