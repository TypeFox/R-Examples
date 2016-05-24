mvProbitLogLikInternal <- function( yMat, xMat, coef, sigma,
   algorithm, nGHK, returnGrad, oneSidedGrad, eps, randomSeed, ... ) {

   # number of regressors
   nReg <- ncol( xMat )

   # checking and preparing model coefficients and correlation coefficients
   coef <- mvProbitPrepareCoef( yMat = yMat, nReg = nReg, coef = coef, 
      sigma = sigma )

   # checking argument 'returnGrad'
   if( length( returnGrad ) != 1 ) {
      stop( "argument 'returnGrad' must be a single logical value" )
   } else if( !is.logical( returnGrad ) ) {
      stop( "argument 'returnGrad' must be logical" )
   }

   # checking argument 'oneSidedGrad'
   if( length( oneSidedGrad ) != 1 ) {
      stop( "argument 'oneSidedGrad' must be a single logical value" )
   } else if( !is.logical( oneSidedGrad ) ) {
      stop( "argument 'oneSidedGrad' must be logical" )
   }

   # checking argument 'eps'
   if( oneSidedGrad ) {
      if( length( eps ) != 1 ) {
         stop( "argument 'eps' must be a single numeric value" )
      } else if( !is.numeric( eps ) ) {
         stop( "argument 'eps' must be numeric" )
      }
   }

   # number of dependent variables
   nDep <- ncol( coef$sigma )

   # number of observations
   nObs <- nrow( xMat )

   # calculating linear predictors
   xBeta <- matrix( NA, nrow = nObs, ncol = nDep )
   for( i in 1:nDep ) {
      xBeta[ , i ] <- xMat %*% coef$betaEq[[ i ]]
   }

   # calculate log likelihood values (for each observation)
   result <- rep( NA, nObs )
   for( i in 1:nObs ){
      ySign <- 2 * yMat[ i, ] - 1
      xBetaTmp <- xBeta[ i, ] * ySign
      sigmaTmp <- diag( ySign ) %*% coef$sigma %*% diag( ySign )
      result[ i ] <- log( pmvnormWrap( upper = xBetaTmp, sigma = sigmaTmp, 
         algorithm = algorithm, nGHK = nGHK, random.seed = randomSeed, ... ) )
   }

   if( returnGrad ) {
      allCoef <- c( coef$beta, coef$sigma[ lower.tri( coef$sigma ) ] )
      grad <- matrix( NA, nrow = length( result ), 
         ncol = length( allCoef ) )
      for( i in 1:nDep ) {
         # gradients of intercepts
         coefLower <- coefUpper <- allCoef
         InterceptNo <- ( i - 1 ) * nReg + 1
         if( oneSidedGrad ) {
            coefUpper[ InterceptNo ] <- allCoef[ InterceptNo ] + eps
            llLower <- result
         } else {
            coefLower[ InterceptNo ] <- allCoef[ InterceptNo ] - eps / 2
            coefUpper[ InterceptNo ] <- allCoef[ InterceptNo ] + eps / 2
            llLower <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
               coef = coefLower, sigma = NULL, 
               algorithm = algorithm, nGHK = nGHK,
               returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
               randomSeed = randomSeed, ... )
         }
         llUpper <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
            coef = coefUpper, sigma = NULL, 
            algorithm = algorithm, nGHK = nGHK,
            returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
            randomSeed = randomSeed, ... )
         grad[ , InterceptNo ] <- ( llUpper - llLower ) / eps
         # gradients of coefficients of other explanatory variables
         if( nReg > 1 ) {
            for( j in 2:nReg ) {
               grad[ , InterceptNo + j - 1 ] <- 
                  grad[ , InterceptNo ] * xMat[ , j ] 
            }
         }
      }
      # gradients of correlation coefficients
      for( i in ( nDep * nReg + 1 ):length( allCoef ) ) {
         coefLower <- coefUpper <- allCoef
         if( oneSidedGrad ) {
            coefUpper[ i ] <- allCoef[ i ] + eps
            llLower <- result
         } else {
            coefLower[ i ] <- allCoef[ i ] - eps / 2
            coefUpper[ i ] <- allCoef[ i ] + eps / 2
            llLower <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
               coef = coefLower, sigma = NULL, 
               algorithm = algorithm, nGHK = nGHK,
               returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
               randomSeed = randomSeed, ... )
         }
         llUpper <- mvProbitLogLikInternal( yMat = yMat, xMat = xMat, 
            coef = coefUpper, sigma = NULL, 
            algorithm = algorithm, nGHK = nGHK,
            returnGrad = FALSE, oneSidedGrad = FALSE, eps = 0, 
            randomSeed = randomSeed, ... )
         grad[ , i ] <- ( llUpper - llLower ) / eps
      }
      colnames( grad ) <- mvProbitCoefNames( nDep = nDep, nReg = nReg )
      attr( result, "gradient" ) <- grad
   }

   return( result )
}
