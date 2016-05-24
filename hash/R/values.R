# ---------------------------------------------------------------------
# values.R
#   values(hash) : returns the values for a hash
# 
# TODO:
#  - Change get to .get in values
#  - na.action
# ---------------------------------------------------------------------

setGeneric( 'values', function(x, ...) standardGeneric( 'values' ) )

setMethod( 'values', 'hash', 
	function(x, keys=NULL, ... ) { 
      if( is.null(keys) ) keys <- keys(x)
      if( ! is.character(keys) ) keys <- make.keys(keys) 
      return(sapply( keys, get, x, ... ))
	}
) 


setGeneric( 'values<-', function(x, ..., value) standardGeneric( 'values<-' ) ) 
setReplaceMethod( 'values', c('hash', 'ANY' ), 
  function(x, ..., value ) {
    keys <- list(...)$keys
    if ( is.null(keys) ) keys <- keys(x) 
    if ( ! is.character(keys) ) keys <- make.keys(keys)

    x[ keys ] <- value
    return(x)
  }
)

# TEST:
# h <- hash( 1:26, letters )
# values(h)
# values(h, keys=1:5 )
# values(h, keys=c(1,1,1:5) )
# values(h, 1:5 )
# values(h, keys=1:5) <- 6:10 
# values(h) <- rev( letters )

