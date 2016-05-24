coef.maxim <- function( object, ... ) {
   return( object$estimate )
}

coef.maxLik <- function( object, ... ) {
   return( object$estimate )
}

coef.summary.maxLik <- function( object, ... ) {
   result <- object$estimate
   return( result )
}
