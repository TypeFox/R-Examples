coef.censReg <- function( object, logSigma = TRUE, ... ) {

   result <- NextMethod( coef, object )
   if( !logSigma ) {
      if( "logSigma" %in% names( result ) ) {
         pos <- length( result )
         if( names( result )[ pos ] != "logSigma" ) {
            stop( "coefficient 'logSigma' must be the last coefficient" )
         }
         result[ pos ] <- exp( result[ pos ] )
         names( result )[ pos ] <- "sigma"
      } else if( "logSigmaMu" %in% names( result ) &&
            "logSigmaNu" %in% names( result ) ) {
         posNu <- length( result )
         posMu <- posNu - 1
         if( names( result )[ posMu ] != "logSigmaMu" ) {
            stop( "coefficient 'logSigmaMu' must be the second-last coefficient" )
         } else if( names( result )[ posNu ] != "logSigmaNu" ) {
            stop( "coefficient 'logSigmaNu' must be the last coefficient" )
         } 
         result[ posMu ] <- exp( result[ posMu ] )
         names( result )[ posMu ] <- "sigmaMu"
         result[ posNu ] <- exp( result[ posNu ] )
         names( result )[ posNu ] <- "sigmaNu"
      } else {
          stop( "the model must contain either coefficient 'logSigma' or",
            " both coefficients 'logSigmaMu' and 'logSigmaNu'" )
      }
   }
   return( result )
}
