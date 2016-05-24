cesRssDeriv <- function( par, yName, xNames, tName, data, vrs, rho1 = NULL,
      rho2 = NULL, rho = NULL, rhoApprox, nested = FALSE, multErr ) {

   # check if coefficients contain NAs
   if( any( is.na( par ) ) ) {
      return( NA )
   }

   # number of exogenous variables
   nExog <- length( xNames )

   withTime <- !is.null( tName )

   # obtain names of coefficients
   coefNames <- cesCoefNames( nExog = nExog, vrs = vrs, 
      returnRho1 = is.null( rho1 ), returnRho2 = is.null( rho2 ), 
      returnRho = is.null( rho ), nested = nested, withTime = withTime )

   # check rhoApprox
   if( !nested ) {
      rhoApprox <- cesCheckRhoApprox( rhoApprox = rhoApprox, withY = TRUE,
         withDeriv = TRUE )
   }

   # add coefficients rho_1, rho_2, and rho, if they are fixed
   par <- cesCoefAddRho( coef = par, vrs = vrs, rho1 = rho1, rho2 = rho2, 
      rho = rho, nExog = nExog, nested = nested )

   # calculate fitted values and residuals
   yHat <- cesCalc( xNames = xNames, tName = tName, data = data, coef = par,
      rhoApprox = rhoApprox[1], nested = nested )
   if( multErr ) {
      resid <- log( data[[ yName ]] ) - log( yHat )
   } else {
      resid <- data[[ yName ]] - yHat
   }

   # obtain derivatives of the CES with respect to coefficients
   derivCoef <- cesDerivCoef( par = par, xNames = xNames, data = data, 
      tName = tName, vrs = vrs, returnRho1 = is.null( rho1 ), 
      returnRho2 = is.null( rho2 ), returnRho = is.null( rho ), 
      rhoApprox = rhoApprox[-1], nested = nested )

   # prepare vector of gradients (to be returned)
   result <- numeric( ncol( derivCoef ) )
   names( result ) <- colnames( derivCoef )
   for( coefName in colnames( derivCoef ) ) {
      if( multErr ) {
         result[ coefName ] <- sum( - 2 * resid * derivCoef[ , coefName ] / yHat )
      } else {
         result[ coefName ] <- sum( - 2 * resid * derivCoef[ , coefName ] )
      }
   }
   return( result )
}
