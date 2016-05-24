mvProbitMargEff <- function( formula, coef, sigma = NULL, vcov = NULL, data,
   cond = FALSE, algorithm = "GHK", nGHK = 1000, eps = 1e-06, 
   dummyVars = NA, addMean = FALSE, returnJacobian = FALSE, 
   random.seed = 123, ... ) {

   # checking argument 'formula'
   if( is.list( formula ) ) {
      stop( "using different regressors for the dependent variables",
         " has not been implemented yet. Sorry!" )
   } else if( class( formula ) != "formula" ) {
      stop( "argument 'formula' must be a formula" )
   }

   # checking argument 'data'
   if( !is.data.frame( data ) ) {
      stop( "argument 'data' must be a data frame" )
   }

   # checking argument 'addMean'
   if( length( addMean ) != 1 ) {
      stop( "argument 'addMean' must be a single logical value" )
   } else if( !is.logical( addMean ) ) {
      stop( "argument 'addMean' must be a logical value" )
   }

   # preparing model matrix
   mc <- match.call( expand.dots = FALSE )
   m <- match( "data", names( mc ), 0 )
   mf <- mc[ c( 1, m ) ]
   mf$formula <- formula
   attributes( mf$formula ) <- NULL
   mf$na.action <- na.pass
   mf[[ 1 ]] <- as.name( "model.frame" )
   mf <- eval( mf, parent.frame() )
   mt <- attr( mf, "terms" )
   xMat <- model.matrix( mt, mf )

   # preparing model response
   yMat <- model.response( mf )
   if( !is.null( yMat ) ) {
      if( !is.matrix( yMat ) ) {
         stop( "either zero or at least two dependent variables",
            " must be specified in argument 'formula'",
            " (e.g. by 'cbind( y1, y2 ) ~ ...')" )
      } else if( !all( yMat %in% c( 0, 1,TRUE, FALSE ) ) ) {
         stop( "all dependent variables must be either 0, 1, TRUE, or FALSE" )
      }
   }

   result <- mvProbitMargEffInternal( yMat = yMat, xMat = xMat, 
      coef = coef, sigma = sigma, cond = cond, 
      algorithm = algorithm, nGHK = nGHK, eps = eps, 
      dummyVars = dummyVars, random.seed = random.seed, ... )

   # join all model coefficients and correlation coefficients
   if( !is.null( sigma ) ) {
      coef <- c( coef, sigma[ lower.tri( sigma ) ] )
   }

   if( !is.null( vcov ) ) {
      # check argument 'vcov'
      if( !is.matrix( vcov ) ) {
         stop( "argument 'vcov' must be a matrix" )
      } else if( nrow( vcov ) != ncol( vcov ) ) {
         stop( "argument 'vcov' must be a quadratic matrix" )
      } else if( !isSymmetric( vcov ) ) {
         stop( "argument 'vcov' must be a symmetric matrix" )
      } else if( nrow( vcov ) != length( coef ) ) {
         stop( "argument 'vcov' must have as many rows and columns",
            " as there are coefficients (model coefficients +",
            " correlation coefficients, i.e. ", length( coef ), ")" )
      }
   }

   # Jacobian matrix d margEff / d coef
   if( !is.null( vcov ) || returnJacobian ) {
      jacobian <- array( NA, 
         c( nrow( result ), ncol( result ), length( coef ) ) )
      for( i in 1:length( coef ) ) {
         coefL <- coefU <- coef
         coefL[ i ] <- coef[ i ] - eps / 2
         coefU[ i ] <- coef[ i ] + eps / 2
         margEffL <- mvProbitMargEffInternal( yMat = yMat, xMat = xMat, 
            coef = coefL, sigma = NULL, cond = cond, 
            algorithm = algorithm, nGHK = nGHK, eps = eps, dummyVars = dummyVars,
            random.seed = random.seed, ... )
         margEffU <- mvProbitMargEffInternal( yMat = yMat, xMat = xMat, 
            coef = coefU, sigma = NULL, cond = cond, 
            algorithm = algorithm, nGHK = nGHK, eps = eps, dummyVars = dummyVars,
            random.seed = random.seed, ... )
         jacobian[ , , i ] <- as.matrix( ( margEffU - margEffL ) / eps )
         coefNames <- names( coef )
         if( is.null( coefNames ) ) {
            coefNames <- mvProbitCoefNames( 
               nDep = round( ncol( result ) / ( ncol( xMat ) - 1 ) ), 
               nReg = ncol( xMat ) )
         }
         dimnames( jacobian ) <- 
            list( rownames( data ), names( result ), coefNames )
      }
   }

   if( !is.null( vcov ) ) {
      margEffCov <- array( NA, 
         c( nrow( result ), ncol( result ), ncol( result ) ) )
      for( i in 1:nrow( result ) ) {
         margEffCov[ i, , ] <- 
            jacobian[ i, , ] %*% vcov %*% t( jacobian[ i, , ] )
      }
      dimnames( margEffCov ) <- 
         list( rownames( data ), names( result ), names( result ) )

      attr( result, "vcov" ) <- margEffCov
   }

   # add mean values of marginal effects if demanded by the user
   if( addMean ) {
      result <- rbind( result, mean = colMeans( result ) )
      if( !is.null( vcov ) || returnJacobian ) {
         mJacobian <- jacobian[ 1, , ]
         if( nrow( xMat ) > 1 ) {
            for( i in 2:nrow( xMat ) ){
               mJacobian <- mJacobian + jacobian[ i, , ]
            }
            mJacobian <- mJacobian / nrow( xMat )
            jacobian <- abind( jacobian, mean = mJacobian, along = 1 )
         }
      }
      if( !is.null( vcov ) ) {
         mVCov <- mJacobian %*% vcov %*% t( mJacobian )
         margEffCov <- abind( margEffCov, mean = mVCov, along = 1 )
      }      
   }

   if( returnJacobian ) {
      attr( result, "jacobian" ) <- jacobian
   }

   if( !is.null( vcov ) ) {
      attr( result, "vcov" ) <- margEffCov
   }

   class( result ) <- c( "mvProbitMargEff", class( result ) )

   return( result )
}
