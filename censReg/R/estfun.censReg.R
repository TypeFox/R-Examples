estfun.censReg <- function( x, includeSigma = TRUE, ... ) {

   result <- NextMethod( estfun, x )

   if( !includeSigma ) {
      result <- result[ , colnames( result ) != "logSigma" ]
   }

   return( result )
}

