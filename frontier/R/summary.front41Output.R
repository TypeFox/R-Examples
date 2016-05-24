summary.front41Output <- function( object, ... ) {

   object$olsResults <- cbind( object$olsResults,
      2 * pt( abs( object$olsResults[ , 3 ] ),
      df = object$nObs - object$nXvars, lower.tail = FALSE ) )
   colnames( object$olsResults )[ 4 ] <- "Pr(>|t|)"

   object$gridResults <- cbind( object$gridResults,
      rep( NA, nrow( object$gridResults ) ) )
   colnames( object$gridResults )[ 4 ] <- "Pr(>|z|)"

   object$mleResults <- cbind( object$mleResults,
      2 * pnorm( abs( object$mleResults[ , 3 ] ), lower.tail = FALSE ) )
   colnames( object$mleResults )[ 4 ] <- "Pr(>|z|)"

   class( object ) <- "summary.front41Output"
   return( object )
}
