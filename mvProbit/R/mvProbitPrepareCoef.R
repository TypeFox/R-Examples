mvProbitPrepareCoef <- function( yMat, nReg, coef, sigma ) {

   # checking argument 'coef'
   if( !is.vector( coef, mode = "numeric" ) ) {
      stop( "argument 'coef' must be a numeric vector" )
   }

   if( !is.null( sigma ) ) {
      # checking argument 'sigma'
      if( !is.matrix( sigma ) ) {
         stop( "argument 'sigma' must be a matrix" )
      } else if( nrow( sigma ) != ncol( sigma ) ) {
         stop( "argument 'sigma' must be a quadratic matrix" )
      } else if( !isSymmetric( sigma ) ) {
         stop( "argument 'sigma' must be a symmetric matrix" )
      } else if( any( abs( diag( sigma ) - 1 ) > 1e-7 ) ) {
         stop( "argument 'sigma' must have ones on its diagonal" )
      } else if( !is.null( yMat ) ) {
         if( ncol( sigma ) != ncol( yMat ) ) {
            stop( "the number of dependent variables specified in argument",
               " 'formula' must be equal to the number of rows and colums",
               " of the matrix specified by argument 'sigma'" )
         }
      }
      # number of dependent variables
      nDep <- ncol( sigma )
      # number of model coefficients
      nCoef <- nDep * nReg
      if( length( coef ) != nCoef ) {
         stop( "given that argument 'sigma' has been specified",
            " argument coef must have ", nCoef, " elements" )
      }
   } else {
      if( !is.null( yMat ) ) {
         # number of dependent variables
         nDep <- ncol( yMat )
         # number of model coefficients
         nCoef <- nDep * nReg
         # number of parameters including sigma
         nCoefSigma <- nCoef + nDep * ( nDep - 1 ) / 2
         if( length( coef ) != nCoefSigma ) {
            stop( "given that argument 'sigma' is 'NULL'",
               " argument coef must have ", nCoefSigma, " elements" )
         }
      } else {
         # number of dependent variables
         nDep <- round( - nReg + 0.5 + 
            sqrt( ( nReg - 0.5 )^2 + 2 * length( coef ) ) )
         # number of model coefficients
         nCoef <- nDep * nReg
         # number of parameters including sigma
         nCoefSigma <- nCoef + nDep * ( nDep - 1 ) / 2
         if( length( coef ) != nCoefSigma ) {
            stop( "given that argument 'sigma' is 'NULL'",
               " argument coef must have ", nCoefSigma, " elements",
               " if the model has ", nDep, " dependent variables" )
         }
      }
      # extracting correlation coefficients from 'coef' if they are there
      sigma <- diag( nDep )
      sigma[ lower.tri( sigma ) ] <- coef[ -c( 1:nCoef ) ]
      sigma[ upper.tri( sigma ) ] <- t( sigma )[ upper.tri( sigma ) ]
      coef <- coef[ 1:nCoef ]
   }

   # separating model coefficients for different equations
   betaEq <- list()
   for( i in 1:nDep ) {
      betaEq[[ i ]] <- coef[ ( ( i - 1 ) * nReg + 1 ):( i * nReg ) ]
   }


   result <- list()
   result$beta <- coef
   result$sigma <- sigma
   result$betaEq <- betaEq

   return( result )
}

