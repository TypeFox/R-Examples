print.summary.TransModel <-
function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("\n")
	printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
	class(x)<-"summary.TransModel"
}
