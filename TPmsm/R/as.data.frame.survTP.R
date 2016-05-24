as.data.frame.survTP <- function(x, ..., package="TPmsm") {
	if ( missing(x) ) stop("Argument 'x' is missing, with no default")
	if ( !is.survTP(x) ) stop("Argument 'x' must be of class 'survTP'")
	package <- match.arg(arg=package, choices=c("TPmsm", "p3state.msm", "etm"), several.ok=FALSE)
	func <- switch(package, "TPmsm"=OutTPmsm, "p3state.msm"=Outp3state, "etm"=Outetm)
	return( func(x[[1]], 1:4) )
}
