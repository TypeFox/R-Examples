addFixedPar <- function( theta, start, fixed, ...) {
   if( is.null( fixed ) ) {
      start <- theta
   } else {
      start[ !fixed ] <- theta
   }

   return( start )
}
