csample.int <- function(n) {
	if ( !is.numeric(n) ) stop("Argument 'n' must be numeric")
	return( .Call(Rf_rsample, as.integer(n), PACKAGE="TPmsm") )
}
