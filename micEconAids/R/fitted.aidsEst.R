fitted.aidsEst <- function( object, ... ) {
   result <- list()
   result$shares <- object$wFitted
   result$quant  <- object$qFitted
   return( result )
}
