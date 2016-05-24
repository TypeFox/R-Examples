# ----------------------------------------------------------------------------------
# METHOD: invert( hash )
#  produces a hash with the values as keys and the keys as values
# ----------------------------------------------------------------------------------

setGeneric( "invert", function(x) standardGeneric( "invert" ) )

setMethod( 'invert', 'hash',
  function(x) {
    h <- hash() 
    for( k in keys(x) ) {
      for( v in make.keys(x[[k]]) ) {
          if ( ! has.key(v,h) ) h[[v]] <- k 
            else h[[v]] <- append( h[[v]], k )
      }
    }

    return(h)
  }

)

# h <- hash( a=1, b=1:2, c=1:3 )
# invert(h)

inverted.hash <- function(...) invert( hash(...) )

# inverted.hash(  a=1, b=1:2, c=1:3 )

