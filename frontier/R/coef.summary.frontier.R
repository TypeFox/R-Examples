coef.summary.frontier <- function( object, which = "mle", ... ) {

   if( tolower( which ) == "ols" ) {
      return( object$olsParam )
   } else if( tolower( which ) == "mle" ) {
      return( object$mleParam )
   } else {
      stop( "argument 'which' must be either 'ols' or 'mle'" )
   }
}