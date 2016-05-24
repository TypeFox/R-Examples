print.HR <- function(x, ...) {
	if ( !inherits(x, "HR") ) stop("Object must be of class HR")
#	mydata <- x$dataset
#	fit <- x$coxfit
#	coxphm <- x$phtest
	cat("\n")
	cat("Cox Proportional Hazards Model\n")
	print(x$coxfit)
	cat("\n")
	cat("Proportional Hazards Assumption\n")
	print(x$phtest)
}
