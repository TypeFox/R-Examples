## Calculate the residual covariance matrix
.calcResidCov <- function( resids, methodResidCov, validObsEq = NULL,
      nCoefEq = NULL, xEq = NULL, diag = FALSE, centered = FALSE,
      useMatrix = FALSE, solvetol = .Machine$double.eps ) {

   eqNames <- NULL
   if( class( resids ) == "data.frame" ) {
      resids <- as.matrix( resids )
      validObsEq <- !is.na( resids )
      eqNames <- names( resids )
   } else if( !is.null( validObsEq ) ) {
      residMat <- matrix( NA, nrow = nrow( validObsEq ),
         ncol = ncol( validObsEq ) )
      for( i in 1:ncol( validObsEq ) ) {
         residMat[ validObsEq[ , i ], i ] <- resids[
            ( 1 + sum( validObsEq[ , 0:(i-1) ] ) ):( sum(validObsEq[ , 1:i ] ) ) ]
      }
      resids <- residMat
      rm( residMat )
   } else {
      stop( "internal error in .calcResidCov: if argument 'validObsEq'",
         " is not provided, argument 'resids' must be a data.frame'" )
   }

   nEq <- ncol( validObsEq )
   result <- matrix( 0, nEq, nEq )
   if( centered ) {
      for( i in 1:nEq ) {
         resids[ , i ] <- resids[ , i ] - mean( resids[ validObsEq[ , i ], i ] )
      }
   }
   validObsAll <- rowSums( !validObsEq ) == 0
   for( i in 1:nEq ) {
      for( j in ifelse( diag, i, 1 ):ifelse( diag, i, nEq ) ) {
         if( methodResidCov == "noDfCor" ) {
            result[ i, j ] <-
               sum( resids[ validObsAll, i ] * resids[ validObsAll, j ] ) /
               sum( validObsAll )
         } else if( methodResidCov == "geomean" ) {
            result[ i, j ] <-
               sum( resids[ validObsAll, i ] * resids[ validObsAll, j ] ) /
               sqrt( ( sum( validObsAll ) - nCoefEq[i] ) * ( sum( validObsAll ) - nCoefEq[j] ) )
         } else if( methodResidCov == "Theil" ) {
            #result[ i, j ] <- sum( residi[[i]] * residi[[j]] ) /
            #   ( sum( validObsAll ) - nCoefEq[i] - nCoefEq[j] + sum( diag(
            #   xEq[[i]] %*% solve( crossprod( xEq[[i]] ), tol=solvetol ) %*%
            #   crossprod( xEq[[i]], xEq[[j]]) %*%
            #   solve( crossprod( xEq[[j]] ), tol=solvetol ) %*%
            #   t( xEq[[j]] ) ) ) )
            result[ i, j ] <-
               sum( resids[ validObsAll, i ] * resids[ validObsAll, j ] ) /
               ( sum( validObsAll ) - nCoefEq[i] - nCoefEq[j] + sum( diag(
               solve( crossprod( xEq[[i]] ), tol=solvetol ) %*%
               crossprod( xEq[[i]], xEq[[j]]) %*%
               solve( crossprod( xEq[[j]] ), tol=solvetol ) %*%
               crossprod( xEq[[j]], xEq[[i]] ) ) ) )

         } else if( methodResidCov == "max" ) {
            result[ i, j ] <-
               sum( resids[ validObsAll, i ] * resids[ validObsAll, j ] ) /
               ( sum( validObsAll ) - max( nCoefEq[ i ], nCoefEq[ j ] ) )
         } else {
            stop( paste( "Argument 'methodResidCov' must be either 'noDfCor',",
                  "'geomean', 'max', or 'Theil'." ) )
         }
      }
   }
   if( !is.null( eqNames ) ) {
      rownames( result ) <- eqNames
      colnames( result ) <- eqNames
   }

   if( useMatrix ){
      result <- as( result, "dspMatrix" )
   }
   return( result )
}

## Calculate Sigma squared
.calcSigma2 <- function( resids, methodResidCov, nObs, nCoef ) {
   if( methodResidCov == "noDfCor" ) {
      result <- sum( resids^2 ) / nObs
   } else if( methodResidCov %in% c( "geomean", "max" ) ){
      result <- sum( resids^2 )/ ( nObs - nCoef )
   } else {
      stop( paste( "Sigma^2 can only be calculated if argument",
         "'methodResidCov' is either 'noDfCor', 'geomean', or 'max'" ) )
   }
}

