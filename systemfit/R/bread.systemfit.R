bread.systemfit <- function ( obj, ... ) {
   if( !is.null( obj$restrict.matrix ) || !is.null( obj$restrict.rhs ) ||
      !is.null( obj$restrict.regMat ) ) {
      stop( "returning the 'bread' for models with restrictions",
            " has not yet been implemented.")
   }
   
   # model matrix
   if( is.null( obj$eq[[1]]$inst ) ) {
      mm <- model.matrix( obj )
   } else {
      mm <- model.matrix( obj, which = "xHat" )
   }
   
   if( is.null( obj$residCovEst ) ) {
      omegaInvXmat <- mm
   } else {
      omegaInvXmat <- t( .calcXtOmegaInv( xMat = mm, sigma = obj$residCovEst, 
         validObsEq = !is.na( residuals( obj ) ), invertSigma = TRUE ) )
   }
   
   result <- solve( crossprod( mm, omegaInvXmat ) / nrow( mm ) )
   
   return( result )
}