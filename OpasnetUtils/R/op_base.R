# Old functions that are now wrappers for the opbase-family for compatibility reasons

op_baseGetData <- function(dsn, ident, ...) {
	return(opbase.data(ident, ...))
}
op_baseGetLocs <- function(dsn, ident, ...) {
	stop('Deprecated method! op_baseGetLocs')
	#opbase.old.locations(dsn, ident, ...)
}