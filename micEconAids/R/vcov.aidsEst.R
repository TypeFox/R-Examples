vcov.aidsEst <- function( object, ... ) {
   result <- object$coef$allcov
   return( result )
}
