print.summary.GORMC <-
function(x,...){
	cat("Call:\n")
	print(x$call)
	cat("\n")
	printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
	print(paste("Loglik=",round(x$loglik,2)))
	class(x)<-"summary.GORMC"
}
