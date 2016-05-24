summary.sparseBC <- function(object, ...){
	cat("Summary for the object \"sparseBC\"\n")
	cat("Call:\n\t"); print(object$cl,...);	cat("\n")


		
	cat("Cluster labels for the rows:\n")
	print(object$Cs)
	cat("\n")
	
	cat("Cluster labels for the columns:\n")
	print(object$Ds)
	cat("\n")
		
	cat("The estimated bicluster means:\n")
	print(object$Mus)
	
	invisible(object)
}
