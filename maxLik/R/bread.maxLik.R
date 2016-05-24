bread.maxLik <- function( x, ... ) {
   return( vcov( x ) * nObs( x ) )
}

