setThreadsTP <-function(num_threads=NULL) {
	if ( is.null(num_threads) ) return( invisible( .Call(Rf_rset_num_threads, num_threads, PACKAGE="TPmsm") ) )
	if ( !is.numeric(num_threads) ) stop("Argument 'num_threads' must be numeric")
	if (length(num_threads) != 1L) stop("Argument 'num_threads' must be of length 1")
	if ( !is.finite(num_threads) ) stop("Argument 'num_threads' must be finite")
	if (num_threads < 1) stop("Argument 'num_threads' must be greater or equal to 1")
	return( invisible( .Call(Rf_rset_num_threads, as.integer(num_threads), PACKAGE="TPmsm") ) )
}
