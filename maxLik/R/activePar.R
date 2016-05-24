## activePar: returns parameters which are free under maximisation (not fixed as constants)

activePar <- function(x, ...)
    UseMethod("activePar")

activePar.default <- function(x, ...) {
   if( !is.null( x$fixed ) ) {
      result <- !x$fixed
   } else {
      result <- x$activePar
   }
   if( is.null( result ) ) {
      result <- rep( TRUE, length( coef( x ) ) )
   }
   return( result )
}
