logLik.systemfit <- function( object, residCovDiag = FALSE, ... ){

   if( length( residCovDiag ) != 1 || any( !is.logical( residCovDiag ) ) ) {
      stop( "argument 'residCovDiag' must be a single logical value" )
   } 
      
   resid <- residuals( object )
   residCov <- .calcResidCov( resid, "noDfCor" )
   if( residCovDiag ) {
      residCov <- diag( diag( residCov ) )
   }
   residCovInv <- solve( residCov )
   resid <- as.matrix( resid )
   nEq <- ncol( resid )

   result <- 0
   for( i in 1:nrow( resid ) ) {
      vEq <- !is.na( resid[ i, ] )
      if( sum( vEq ) == nEq ) {
         result <- result - ( nEq / 2 ) * log( 2 * pi ) -
            ( 1 / 2 ) * log( det( residCov ) ) -
            ( 1 / 2 ) * resid[ i, , drop = FALSE ] %*% residCovInv %*%
               t( resid[ i, , drop = FALSE ] )
      } else if( sum( vEq ) > 0 ) {
         nEq2 <- sum( vEq )
         residCov2 <- residCov[ vEq, vEq, drop = FALSE ] -
            residCov[ vEq, !vEq, drop = FALSE ] %*%
            solve( residCov[ !vEq, !vEq, drop = FALSE ] ) %*%
            residCov[ !vEq, vEq, drop = FALSE ]
         residCov2Inv <- solve( residCov2 )
         result <- result - ( nEq2 / 2 ) * log( 2 * pi ) -
            ( 1 / 2 ) * log( det( residCov2 ) ) -
            ( 1 / 2 ) * resid[ i, vEq, drop = FALSE ] %*%
               residCov2Inv %*% t( resid[ i, vEq, drop = FALSE ] )
      }
   }

   if( object$method %in% c( "OLS", "2SLS" ) ){
      nSigma <- 1
   } else if( object$method %in% c( "WLS", "W2SLS" ) ){
      nSigma <- nEq
   } else if( object$method %in% c( "SUR", "3SLS" ) ){
      nSigma <- nEq * ( nEq + 1 ) / 2
   } else {
      stop( "internal error: unknown estimation method '", object$method, "'" )
   }

   attributes( result )$nobs <- df.residual( object ) + object$rank
   attributes( result )$df <- object$rank + nSigma
   class( result ) <- "logLik"

   return( result )
}