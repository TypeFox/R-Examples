#' @rdname pcps.curve
#' @encoding UTF-8
#' @export
print.pcpscurve<-function(x, ...){
	res<-summary(x)
	cat("Call:\n")
	cat(deparse(x$call), "\n\n")
	cat("PCPS curve observed:\n")
	print(res$curve.obs)
	if(!is.null(x$curve.null.ts)){
		cat("\n")
		cat("Mean PCPS curve null TS:\n")
		print(as.matrix(res$null.model.ts[,1:2]))
	}
	if(!is.null(x$curve.null.bm)){
		cat("\n")
		cat("Mean PCPS curve null BM:\n")
		print(as.matrix(res$null.model.bm[,1:2]))
	}
	invisible(x)
}