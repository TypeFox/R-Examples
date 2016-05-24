mvProbitMargEffInternal <- function( yMat, xMat, coef, sigma,
   cond, algorithm, nGHK, eps, dummyVars,
   random.seed, ... ) {

   # number of observations
   nObs <- nrow( xMat )

   # number of regressors
   nReg <- ncol( xMat )

   # names of regressors
   xNames <- colnames( xMat )

   # detect dummy variables if they should be determined automatically
   if( !is.null( dummyVars ) && is.na( dummyVars[ 1 ] ) ) {
      dummyVars <- NULL
      for( i in 2:nReg ) {
         if( all( xMat[ , i ] %in% c( 0, 1, FALSE, TRUE ) ) ) {
            dummyVars <- c( dummyVars, xNames[ i ] )
         }
      }
   }

   # check dummy variables
   if( !is.null( dummyVars ) ) {
      foundDummyVars <- dummyVars %in% c( "", xNames[ - 1 ] )
      if( ! all( foundDummyVars ) ) {
         warning( "variable(s) '", 
            paste( dummyVars[ !foundDummyVars ], collapse = "', '" ),
            "' specified in argument 'dummyVars'",
            " seem(s) to be no explanatory variable(s)" )
      }
   }

   # checking and preparing model coefficients and correlation coefficients
   coef <- mvProbitPrepareCoef( yMat = yMat, nReg = nReg, coef = coef, 
      sigma = sigma )

   # number of dependent variables
   nDep <- ncol( coef$sigma )

   # calculate marginal effects of X'beta
   margEffXBeta <- array( NA, c( nObs, nDep, nDep ) )
   if( !all( xNames %in% dummyVars ) ) {
      for( i in 1:nDep ) {
         betaLower <- betaUpper <- coef$beta
         interceptNo <- ( i - 1 ) * nReg + 1
         betaLower[ interceptNo ] <- coef$beta[ interceptNo ] - eps / 2
         betaUpper[ interceptNo ] <- coef$beta[ interceptNo ] + eps / 2
         yMatL <- mvProbitExpInternal( yMat = yMat, xMat = xMat,
            coef = betaLower, sigma = coef$sigma, 
            cond = cond, algorithm = algorithm, nGHK = nGHK, 
            random.seed = random.seed, ... )
         yMatU <- mvProbitExpInternal( yMat = yMat, xMat = xMat, 
            coef = betaUpper, sigma = coef$sigma, 
            cond = cond, algorithm = algorithm, nGHK = nGHK, 
            random.seed = random.seed, ... )
         margEffXBeta[ , , i ] <- as.matrix( yMatU - yMatL ) / eps
      }
   }

   # calculate marginal effects
   for( i in 2:nReg ) {
      xMatL <- xMatU <- xMat
      isDummy <- xNames[ i ] %in% dummyVars
      if( isDummy ) {
         xMatL[ , i ] <- 0
         xMatU[ , i ] <- 1
         yMatL <- mvProbitExpInternal( yMat = yMat, xMat = xMatL,
            coef = coef$beta, sigma = coef$sigma, 
            cond = cond, algorithm = algorithm, nGHK = nGHK, 
            random.seed = random.seed, ... )
         yMatU <- mvProbitExpInternal( yMat = yMat, xMat = xMatU, 
            coef = coef$beta, sigma = coef$sigma, 
            cond = cond, algorithm = algorithm, nGHK = nGHK, 
            random.seed = random.seed, ... )
         margEff <- yMatU - yMatL
      } else {
         margEff <- matrix( 0, nrow = nObs, ncol = nDep )
         for( j in 1:nDep ) {
            coefNo <- ( j - 1 ) * nReg + i
            margEff <- margEff + coef$beta[ coefNo ] * margEffXBeta[ , , j ]
         }
         margEff <- as.data.frame( margEff )
      }

      # names of dependent variables
      if( i == 2 ) {
         if( !is.null( yMat ) ) {
            yNames <- colnames( yMat )
         } else {
            yNames <- paste( "y", 1:ncol( margEff ), sep = "" )
         }
      }

      # label marginal effects
      names( margEff ) <- paste( "d", yNames, "d", xNames[ i ], sep = "_" )

      if( i == 2 ) {
         result <- margEff
      } else {
         result <- cbind( result, margEff )
      }
   }

   return( result )
}
