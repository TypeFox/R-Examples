summary.analytics <- function(object, ...) {	
	algorithms <- object@algorithm_summary
	average <- colSums(algorithms, na.rm=T)/nrow(algorithms)
	
	cat("ENSEMBLE SUMMARY\n\n")
	print(object@ensemble_summary)
	cat("\n\n")
	cat("ALGORITHM PERFORMANCE\n\n")
	print(average)
}