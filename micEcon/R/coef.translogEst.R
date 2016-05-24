coef.translogEst <- function( object, ... ) {
   return( coef.quadFuncEst( object ) )
}

vcov.translogEst <- function( object, ... ) {
   return( vcov.quadFuncEst( object ) )
}
