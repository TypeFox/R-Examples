mvProbitExpInternal <- function( yMat, xMat, coef, sigma,
   cond, algorithm, nGHK, random.seed, ... ) {

   # checking argument 'yMat'
   if( !is.null( yMat ) ) {
      if( !all( yMat %in% c( 0, 1,TRUE, FALSE ) ) ) {
         stop( "all dependent variables must be either 0, 1, TRUE, or FALSE" )
      }
   }

   # checking argument 'cond'
   if( !is.logical( cond ) ) {
      stop( "argument 'cond' must be logical" )
   } else if( length( cond ) != 1 ) {
      stop( "argument 'cond' must be a single logical values" )
   }

   # number of regressors
   nReg <- ncol( xMat )

   # number of observations
   nObs <- nrow( xMat )

   # checking and preparing model coefficients and correlation coefficients
   coef <- mvProbitPrepareCoef( yMat = yMat, nReg = nReg, coef = coef, 
      sigma = sigma )

   # number of dependent variables
   nDep <- ncol( coef$sigma )

   # calculating linear predictors
   xBeta <- matrix( NA, nrow = nObs, ncol = nDep )
   for( i in 1:nDep ) {
      xBeta[ , i ] <- xMat %*% coef$betaEq[[ i ]]
   }

   if( cond ) {
      # conditional expectations
      result <- matrix( NA, nrow = nObs, ncol = nDep )
      if( is.null( yMat ) ) {
         # assuming that all other dependent variables are one
         for( i in 1:nObs ) {
            for( k in 1:nDep ) {
               result[ i, k ] <- 
                  pmvnormWrap( upper = xBeta[ i, ], sigma = coef$sigma,
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... ) /
                  pmvnormWrap( upper = xBeta[ i, -k ],
                     sigma = coef$sigma[ -k, -k, drop = FALSE ],
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... )
            }
         }
      } else {
         # assuming that all other dependent variables are as observed
         for( i in 1:nObs ){
            for( k in 1:nDep ) {
               ySign <- 2 * yMat[ i, ] - 1
               ySign[ k ] <- 1
               xBetaTmp <- xBeta[ i, ] * ySign
               sigmaTmp <- diag( ySign ) %*% coef$sigma %*% diag( ySign )
               result[ i, k ] <- 
                  pmvnormWrap( upper = xBetaTmp, sigma = sigmaTmp,
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... ) / 
                  pmvnormWrap( upper = xBetaTmp[ -k ],
                     sigma = sigmaTmp[ -k, -k, drop = FALSE ],
                     algorithm = algorithm, nGHK = nGHK, 
                     random.seed = random.seed, ... )
            }
         }
      }
   } else {
      result <- pnorm( xBeta )
   }

   if( !is.null( yMat ) ) {
      colnames( result ) <- colnames( yMat )
   }

   result <- as.data.frame( result )

   return( result )
}
