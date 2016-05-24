summary.edr <-
function(object,...){
	x <- object
	if (!inherits(x, "edr")) 
		stop("use only with \"edr\" objects")

	cat(paste("Reduction method performed:", x$method),"\n")
	cat(" \n")

	cat(paste("Number of observations:", x$n),"\n")
	cat(paste("Dimension reduction K:", x$K),"\n")
	cat(paste("Number of slices:", paste(x$H, collapse=", ")),"\n")
	cat(" \n")

	if (is.null(x$matEDR)) {
		cat("Indices estimation results:\n")
		tmp <- x$indices
		colnames(tmp) <- 1:x$K
		row.names(tmp) <- paste("estimated index", 1:x$n, sep=" ") 
	} else {
		cat("Result of EDR directions estimation:\n" )
		tmp <- matrix(x$matEDR[,1:x$K],ncol=x$K)
		colnames(tmp) <- paste("estimated direction",1:x$K,sep=" ")
		row.names(tmp) <- 1:dim(tmp)[1]
	}	
	prmatrix(signif(tmp,3))
	cat("\n")

	if (!is.null(x$eigvalEDR)) {
		cat("List of eigenvalues: \n" )	
		prmatrix(signif(x$eigvalEDR,3))
		cat("\n")
}	}

