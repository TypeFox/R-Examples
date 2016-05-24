summary.cesEst <- function( object, rSquaredLog = object$multErr,
      ela = TRUE, ... ) {

   # number of observations
   nObs <- length( residuals( object ) )

   # square root of the estimated variance of the random error
   object$sigma <- sqrt( object$rss / nObs )

   # R-squared value
   fitVal <- fitted( object ) # the returned fitted values are never logged
   if( rSquaredLog ) {
      if( object$multErr ) {
         yVal <- log( fitVal ) + residuals( object )  # log( y )
      } else {
         yVal <- log( fitVal + residuals( object ) )  # log( y )
      }
      res <- yVal - log( fitVal )
   } else {
      if( object$multErr ) {
         yVal <- fitVal * exp( residuals( object ) )  # y
      } else {
         yVal <- fitVal + residuals( object )         # y
      }
      res <- yVal - fitVal
   }
   object$r.squared <- rSquared( y = yVal, resid = res )

   # covariance matrix of the estimated coefficients/parameters
   if( is.null( object$vcov ) ) {
      object$vcov <- object$sigma^2 * object$cov.unscaled
   }

   if( ela && !is.null( object$ela ) ) {
      # (co)variances of the elasticities of substitution
      elaGrad <- matrix( 0, nrow = length( object$ela ), 
         ncol = length( coef( object ) ),
         dimnames = list( names( object$ela ), names( coef( object ) ) ) )
      elaGrad[ nrow( elaGrad ), "rho" ] <- - 1 / ( 1 + coef( object )[ "rho" ] )^2
      if( nrow( elaGrad ) >= 2 ) {
         elaGrad[ 1, "rho_1" ] <- - 1 / ( 1 + coef( object )[ "rho_1" ] )^2
      }
      if( nrow( elaGrad ) == 3 ) {
         elaGrad[ 2, "rho_2" ] <- - 1 / ( 1 + coef( object )[ "rho_2" ] )^2
      }
      object$elaCov <- elaGrad %*% object$vcov %*% t( elaGrad )
      object$elaCov[ is.na( object$ela ), ] <- NA
      object$elaCov[ , is.na( object$ela ) ] <- NA
   }

   object$coefficients <- coefTable( coef( object ),
      diag( object$vcov )^0.5, df = Inf )

   if( ela && !is.null( object$elaCov ) ) {
      object$ela <- coefTable( object$ela, diag( object$elaCov )^0.5, df = Inf )
   } else {
      object$ela <- NULL
   }


   class( object ) <- c( "summary.cesEst", class( object ) )
   return( object )
}
