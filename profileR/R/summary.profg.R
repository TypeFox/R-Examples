
#'@method summary profg
#'@export

summary.profg <- function(object, ...){
	cat("Call:\n")
	print(object$call)
	cat("\nHypothesis Tests:\n")
	print(object$profile.test)
}

