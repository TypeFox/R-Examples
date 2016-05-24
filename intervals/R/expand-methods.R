# Not exported
adjust <- function( x, delta, type, direction = 1 ) {
  signs <- rep( c( direction, -direction ), c( nrow(x), nrow(x) ) )
  if ( nrow(x) %% length(delta) != 0 )
    warning( "Length of 'delta' does not evenly divide number of intervals.", call. = FALSE )
  delta <- rep( delta, length = nrow( x ) )
  switch(
         type,
         relative = x * ( 1 - delta ) ^ ( sign(x) * signs ),
         absolute = x - delta * signs
         )
}

setGeneric( "expand", def = function( x, ... ) standardGeneric( "expand" ) )

setMethod(
          "expand",
          signature( "Intervals_virtual" ),
          function(
                   x,
                   delta = 0,
                   type = c( "absolute", "relative" )
                   )
          {
            if ( any( delta < 0, na.rm = TRUE ) )
              stop( "The 'delta' argument should not contain negative values." )
            type = match.arg( type )
            if ( type(x) == "Z" && ( type == "relative" || any( delta %% 1 != 0, na.rm = TRUE ) ) )
              stop( "Only absolute, integer-valued expansion permitted for type 'Z'." )
            x@.Data <- adjust( x, delta, type, 1 )
            return( x )
          }
          )

setGeneric( "contract", def = function( x, ... ) standardGeneric( "contract" ) )

setMethod(
          "contract",
          signature( "Intervals_virtual" ),
          function(
                   x,
                   delta = 0,
                   type = c( "absolute", "relative" )
                   )
          {
            if ( any( delta < 0, na.rm = TRUE ) )
              stop( "The 'delta' argument should not contain negative values." )
            type = match.arg( type )
            if ( type(x) == "Z" && ( type == "relative" || any( delta %% 1 != 0, na.rm = TRUE ) ) )
              stop( "Only absolute, integer-valued contraction permitted for type 'Z'." )
            x@.Data <- adjust( x, delta, type, -1 )
            drop <- x[,1] > x[,2] | empty(x)
            if ( any( drop ) ) {
              warning( "Some empty intervals eliminated.", call. = FALSE )
              x <- x[ !drop, ]
            }
            return( x )
          }
          )
