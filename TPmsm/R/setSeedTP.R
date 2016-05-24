setSeedTP <- function(x) {
	UseMethod("setSeedTP")
}

setSeedTP.NULL <- function(x) {
	if ( missing(x) ) x <- NULL
	if ( !is.null(x) ) stop("Argument 'x' must be NULL")
	return( invisible( .Call(Rf_rset_seed, x, PACKAGE="TPmsm") ) )
}

setSeedTP.TPmsmSeed <- function(x) {
	if ( !is.TPmsmSeed(x) ) stop("Argument 'x' must be of class 'TPmsmSeed'")
	if ( length(x) != setThreadsTP(NULL) ) stop("Object's length must match predefined number of threads")
	return( invisible( .Call(Rf_rset_seed, x, PACKAGE="TPmsm") ) )
}

is.TPmsmSeed <- function(x) {
	ret <- inherits(x, "TPmsmSeed") & is.list(x) & (length(x) >= 1)
	if (!ret) return(ret)
	for ( i in 1:length(x) ) {
		ret <- ret & is.double(x[[i]])
		ret <- ret & (length(x[[i]]) == 6)
	}
	return(ret)
}
