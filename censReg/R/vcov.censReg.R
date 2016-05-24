vcov.censReg <- function( object, logSigma = TRUE, ... ) {

   result <- NextMethod( vcov, object )
   if( !logSigma ) {
      if( "logSigma" %in% rownames( result ) ) {
         pos <- nrow( result )
         if( rownames( result )[ pos ] != "logSigma" ) {
            stop( "coefficient 'logSigma' must be the last coefficient" )
         }
         jac <- diag( pos )
         jac[ pos, pos ] <- coef( object, logSigma = FALSE )[ pos ]
         dNames <- dimnames( result )
         result <- t( jac ) %*% result %*% jac
         dimnames( result ) <- dNames
         rownames( result )[ pos ] <- colnames( result )[ pos ] <- "sigma"
      } else if( "logSigmaMu" %in% rownames( result ) &&
            "logSigmaNu" %in% rownames( result ) ) {
         posNu <- nrow( result )
         posMu <- posNu - 1
         if( rownames( result )[ posMu ] != "logSigmaMu" ) {
            stop( "coefficient 'logSigmaMu' must be the second-last coefficient" )
         } else if( rownames( result )[ posNu ] != "logSigmaNu" ) {
            stop( "coefficient 'logSigmaNu' must be the last coefficient" )
         }
         jac <- diag( posNu )
         jac[ posMu, posMu ] <- coef( object, logSigma = FALSE )[ posMu ]
         jac[ posNu, posNu ] <- coef( object, logSigma = FALSE )[ posNu ]
         dNames <- dimnames( result )
         result <- t( jac ) %*% result %*% jac
         dimnames( result ) <- dNames
         rownames( result )[ posMu ] <- colnames( result )[ posMu ] <- "sigmaMu"
         rownames( result )[ posNu ] <- colnames( result )[ posNu ] <- "sigmaNu"
      } else {
          stop( "the model must contain either coefficient 'logSigma' or",
            " both coefficients 'logSigmaMu' and 'logSigmaNu'" )
      }
   }
   return( result )
}
