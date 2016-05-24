coef.front41Output <- function( object, which = "MLE", ... ) {

   if( which %in% c( "OLS", "ols" ) ) {
      result <- drop( object$olsResults[ , 1 ] )
   } else if( which %in% c( "GRID", "Grid", "grid" ) ) {
      result <- drop( object$gridResults[ , 1 ] )
   } else if( which %in% c( "MLE", "mle" ) ) {
      result <- drop( object$mleResults[ , 1 ] )
   } else {
      stop( "argument 'which' must be either 'OLS', 'GRID', or 'MLE'" )
   }

   return( result )
}

coef.summary.front41Output <- function( object, which = "MLE", ... ) {

   if( which %in% c( "OLS", "ols" ) ) {
      result <- object$olsResults
   } else if( which %in% c( "GRID", "Grid", "grid" ) ) {
      result <- object$gridResults
   } else if( which %in% c( "MLE", "mle" ) ) {
      result <- object$mleResults
   } else {
      stop( "argument 'which' must be either 'OLS', 'GRID', or 'MLE'" )
   }

   return( result )
}

vcov.front41Output <- function( object, ... ) {
   return( object$mleCov )
}
