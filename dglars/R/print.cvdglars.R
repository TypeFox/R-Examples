print.cvdglars <- function (x,digits = max(3, getOption("digits") - 3), ...){
	b <- x$beta[abs(x$beta)>0]
	cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
	print.default(format(b, digits = digits),print.gap = 2, quote = FALSE,row.names=FALSE)
	cat("\nCross-validation deviance = ",format(min(x$dev_m), digits = digits),"( n. fold = ",x$control$nfold,")")
	cat("\nAlgorithm", x$control$algorithm,"( method =",x$control$method,") with exit =",x$conv,"\n\n")
	invisible(b)
}
