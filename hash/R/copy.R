# -----------------------------------------------------------------------------
# copy.R
#
# -----------------------------------------------------------------------------
setGeneric( 'copy', function(x,...) standardGeneric( 'copy'  ) )

setMethod( 'copy', 'hash', 
  function(x, ... ) {
    
    hash( mget( keys(x), x@.xData ) ) 

  }
)



