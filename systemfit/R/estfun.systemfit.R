estfun.systemfit <- function ( obj, residFit = TRUE, ... ) {
   if( !is.null( obj$restrict.matrix ) || !is.null( obj$restrict.rhs ) ||
        !is.null( obj$restrict.regMat ) ) {
      stop( "returning the estimation function for models with restrictions",
            " has not yet been implemented.")
   }
   
   # residuals
   res <- unlist(  residuals( obj ) )
   
   # model matrix
   if( is.null( obj$eq[[1]]$inst ) ) {
      mm <- model.matrix( obj )
   } else {
      mm <- model.matrix( obj, which = "xHat" )
      if( residFit ) {
         res[ !is.na( res ) ] <- res[ !is.na( res ) ] +
            ( model.matrix( obj ) - mm ) %*% coef( obj )
         #       resid_fit = y - x_fit b
         #       resid = y - x b
         #       resid_fit - resid = - x_fit b + x b
         #       resid_fit = resid + ( x - x_fit ) b
      }
   }
   
   if( sum( !is.na( res ) ) != nrow( mm ) ) {
      stop( "internal error: the number of residuals is not equal to the",
         " number of rows of the model matrix. Please contact the maintainer." )
   }

   if( is.null( obj$residCovEst ) ) {
      omegaInvXmat <- mm
   } else {
      omegaInvXmat <- t( .calcXtOmegaInv( xMat = mm, sigma = obj$residCovEst, 
         validObsEq = !is.na( residuals( obj ) ), invertSigma = TRUE ) )
   }
   
   result <- res[ !is.na( res ) ] * omegaInvXmat
   
   dimnames( result ) <- dimnames( mm )
   
   if( max( abs( colSums( result ) ) ) > 1e-6 ) {
      warning( "the columns of the returned estimating function",
         " do not all sum up to zero,",
         " which indicates that the wrong estimating function is returned" )
   }
   
   return( result )
}