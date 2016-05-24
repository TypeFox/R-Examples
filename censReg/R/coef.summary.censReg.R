coef.summary.censReg <- function( object, logSigma = TRUE, ... ) {

   result <- NextMethod( coef, object )
   if( !logSigma ) {
      if( "logSigma" %in% rownames( result ) ) {
         pos <- nrow( result )
         if( rownames( result )[ pos ] != "logSigma" ) {
            stop( "coefficient 'logSigma' must be the last coefficient" )
         }
         result[ pos, 1 ] <- exp( result[ pos, 1 ] )
         result[ pos, 2 ] <- result[ pos, 1 ] * result[ pos, 2 ]
         result[ pos, 3 ] <- result[ pos, 1 ] / result[ pos, 2 ]
         result[ pos, 4 ] <- 2 * pnorm( abs( result[ pos, 3 ] ), 
            lower.tail = FALSE )
         rownames( result )[ pos ] <- "sigma"
      } else if( "logSigmaMu" %in% rownames( result ) &&
            "logSigmaNu" %in% rownames( result ) ) {
         posNu <- nrow( result )
         posMu <- posNu - 1
         if( rownames( result )[ posMu ] != "logSigmaMu" ) {
            stop( "coefficient 'logSigmaMu' must be the second-last coefficient" )
         } else if( rownames( result )[ posNu ] != "logSigmaNu" ) {
            stop( "coefficient 'logSigmaNu' must be the last coefficient" )
         } 
         result[ posMu:posNu, 1 ] <- exp( result[ posMu:posNu, 1 ] )
         result[ posMu:posNu, 2 ] <- result[ posMu:posNu, 1 ] * 
            result[ posMu:posNu, 2 ]
         result[ posMu:posNu, 3 ] <- result[ posMu:posNu, 1 ] / 
            result[ posMu:posNu, 2 ]
         result[ posMu:posNu, 4 ] <- 2 * 
            pnorm( abs( result[ posMu:posNu, 3 ] ), lower.tail = FALSE )
         rownames( result )[ posMu:posNu ] <- c( "sigmaMu", "sigmaNu" )
      } else {
          stop( "the model must contain either coefficient 'logSigma' or",
            " both coefficients 'logSigmaMu' and 'logSigmaNu'" )
      }
   }
   return( result )
}
