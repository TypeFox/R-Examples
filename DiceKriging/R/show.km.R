show.km <- function(object) {	
	cat("\n")
	cat("Call:\n")
	print(object@call)
	
	p <- length(object@trend.coef)
	trend.coef <- matrix(object@trend.coef, p, 1)
	trend.coef <- t(formatC(trend.coef, width=10, digits=4, format="f", flag=" "))
	trend.names <- formatC(colnames(object@F), width=12)
	
	if ((identical(object@known.param, "Trend")) | (identical(object@known.param, "All"))) {
		dimnames(trend.coef) <- list("", trend.names)
	} else {
		dimnames(trend.coef) <- list("  Estimate", trend.names)
	}
	cat("\n")
	cat("Trend  coeff.:\n")
	print(t(trend.coef), quote=FALSE) 
	
	show(object@covariance)

}