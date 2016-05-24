.ftest.systemfit <- function( object, restrict.matrix,
   restrict.rhs = NULL, vcov. = NULL ){

   coef <- coef( object )

   # coefficient covariance matrix
   if( is.null( vcov. ) ){
      vcov <- vcov( object )
   } else if( is.function( vcov. ) ){
      vcov <- vcov.( object )
   } else {
      vcov <- vcov.
   }

   resid <- unlist( residuals( object ) )
   resid <- resid[ !is.na( resid ) ]
   nEq   <- length( object$eq )
   if( is.null( object$residCovEst ) ) {
      rcov <- diag( nEq )
   } else {
      rcov <- object$residCovEst
   }

   validObsEq <- matrix( NA, nrow = nrow( residuals( object ) ), ncol = nEq )
   for( i in 1:nEq ) {
      validObsEq[ , i ] <- !is.na( residuals( object$eq[[ i ]] ) )
   }

   result <- list()

   result$nRestr <- nrow( restrict.matrix )
   result$df.residual.sys  <- object$df.residual

   numerator <- t( restrict.matrix %*% coef - restrict.rhs ) %*%
      solve( restrict.matrix %*% vcov %*% t( restrict.matrix ) ) %*%
      ( restrict.matrix %*% coef - restrict.rhs )

   resid <- matrix( resid, ncol = 1 )
   if( object$control$useMatrix ) {
      resid <- as( resid, "dgCMatrix" )
      rcov <- as( rcov, "dspMatrix" )
   }

   denominator <- as.numeric( .calcXtOmegaInv( xMat = resid, sigma = rcov,
      validObsEq = validObsEq, useMatrix = object$control$useMatrix ) %*% resid )

   result$statistic <- ( numerator / result$nRestr ) /
      ( denominator / result$df.residual.sys )

   result$p.value <- 1 - pf( result$statistic, result$nRestr, result$df.residual.sys )

   return( result )
}
