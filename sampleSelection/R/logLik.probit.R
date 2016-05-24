logLik.probit <- function( object, ... ) {
   if( !inherits( object, "maxLik" ) ) {
      stop( "the object must inherit from class 'maxLik'" )
   }
   obj <- object
   class( obj ) <- "maxLik"
   result <- logLik( obj )
   attr( result, "nobs" ) <- nObs( object )
   attr( result, "df" ) <- sum( activePar( object ) )
   class( result ) <- "logLik"
   return( result )
}
