#' @rdname pcoa.sig
#' @encoding UTF-8
#' @export
print.pcoasig<-function(x , ...){
	cat("Call:\n")
	cat(deparse(x$call), "\n\n")
	cat("PCoA values:\n")
	print(as.matrix(x$PCoA$values))
	cat("\nProbabilities:\n")
	print(as.matrix(x$probabilities))
	invisible(x)
}