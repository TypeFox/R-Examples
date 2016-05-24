# -----------------------------------------------------------------------------
# keys.R
# METHOD: keys
# -----------------------------------------------------------------------------
setGeneric( "keys", function(x) standardGeneric("keys") )
setMethod( "keys" , "hash" ,
	function(x) ls(x@.Data, all.names=T )
)

names.hash <- function(x) keys(x)


