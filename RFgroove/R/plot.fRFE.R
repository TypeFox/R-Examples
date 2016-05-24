plot.fRFE <-
function(x, ...){
	if(class(x)!="fRFE")
		stop("Wrong class")

	if(is.null(x$histSel[[1]]$indexIter)){
		plot(1:length(x$error), x$error, type="b", xlab="Number of variables", ylab="Validation error estimate", ...)
		abline(v=x$nselected, col="grey", lty=2)
	}else{
		plot(x$histSel[[1]]$indexIter, x$error, type="b", xlab="Number of variables", ylab="Validation error estimate", ...)
		abline(v=x$nselected, col="grey", lty=2)
	}
}
