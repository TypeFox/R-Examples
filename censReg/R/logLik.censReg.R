logLik.censReg <- function( object, ... ) {

   result <- NextMethod( logLik, object )
   attr( result, "df" ) <- sum( activePar( object ) )

   class( result ) <- "logLik"

   return( result )
}
