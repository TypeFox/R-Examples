summary.fRFE <-
function(object, ...){
	if(class(object)!="fRFE")
		stop("Wrong class")

	cat("\n\n --- \t\t\t\t ---")
	cat("\n --- \tSummary functional RFE\t ---")
	cat("\n --- \t\t\t\t ---\n\n\n")

	if(!is.null(object$nclust) & !is.null(object$K))
		cat("\nStability selection using", object$K, "subsamples, parallel execution on", object$nclust, "cores\n\n")

	cat("Number of selected variables using a validation set:", object$nselected, "\n\n")

	cat("Selected variables:\n")
	cat(object$selection, "\n\n")

	cat("Validation error for the best model:", round(min(object$error), 4), "\n\n")
}
