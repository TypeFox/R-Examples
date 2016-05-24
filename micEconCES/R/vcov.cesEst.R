vcov.cesEst <- function( object, ... ){

   if( !is.null( object$vcov ) ) {
      return( object$vcov )
   } else {
      return( summary( object )$vcov )
   }
}

