#' Hash/Dictionary Lookup
#' 
#' \code{\%l*\%} - A deprecated binary operator version of \code{lookup}.  This
#' will be removed in a subsequent version of \pkg{qdapTools}.  Use \code{\%l\%}
#' instead.
#' 
#' @export
#' @rdname Deprecated
`%l*%` <- function(terms, key.match) {
    .Deprecated(msg = paste("`%l*%` is deprecated.  Please use `%l%` instead."), 
    	old = as.character(sys.call(sys.parent()))[1L])
	
    lookup(terms = terms, key.match = key.match, missing = NULL)
}