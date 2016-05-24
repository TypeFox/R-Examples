snqProfitImposeConvexity <- function( estResult, rankReduction = 0,
   start = 10, optimMethod = "BFGS", control = list( maxit=5000 ),
   stErMethod = "none", nRep = 1000, verbose = 0 ) {

   if( class( estResult ) != "snqProfitEst" ) {
      stop( "argument 'estResult' must be of class 'snqProfitEst'" )
   }
   if( !( stErMethod %in% c( "none", "jackknife", "resample", "coefSim" ) ) ) {
      stop( "argument 'stErMethod' must be either 'none', 'resample',",
         " 'jackknife' or 'coefSim'" )
   }
   if( estResult$convexity ) {
      warning( "This profit function is already convex in prices" )
      return( estResult )
   }
   priceNames  <- estResult$priceNames
   quantNames  <- estResult$quantNames
   fixNames  <- estResult$fixNames
   nNetput <- length( priceNames )
   nFix    <- length( fixNames )
   nObs    <- nrow( estResult$data )
   result  <- list()

   ## preparations for minimum distance estimation
   uVecliHessian <- vecli( estResult$hessian[ 1:( nNetput - 1 ),
      1:( nNetput - 1 ) ] )
      # vector of linear indep. values of unconstr. Hessian
   hessianDeriv <- snqProfitHessianDeriv( estResult$pMeans, estResult$weights,
      nFix = nFix, form = estResult$form )
      # derivatives of the Hessian with respect to the coefficients

   ## distance function between the unconstrained and constrained Hessian
   hessianDistance <- function( cholVec ) {
         # cholVec = vector of values of triangular Cholesky Matrix
      cholMat <- triang( cholVec, nNetput - 1 )     # triangular Cholesky Matrix
      cVecliHessian <- vecli( t(cholMat) %*% cholMat )
         # vector of linear independent values of constrained Hessian
      result <- t( uVecliHessian - cVecliHessian ) %*% solve( hessianDeriv %*%
         estResult$coef$liCoefCov %*% t( hessianDeriv ) ) %*%
         ( uVecliHessian - cVecliHessian )
      return( result )
   }

   ## non-linear minimization of the distance function
   if( verbose >= 1 ) {
      cat( "Starting minimization of the distance function\n" )
   }
   nCholValues <- nNetput * ( nNetput - 1 ) / 2 -
      rankReduction * ( rankReduction + 1 ) / 2
      # number of non-zero values in the Cholesky matrix
   cholVec <- array( start, c( nCholValues ) ) # starting values
   mindist <- optim( cholVec, hessianDistance, method = optimMethod,
      control = control )
   if( mindist$convergence != 0 ) {
      if( mindist$convergence == 1 ) {
         stop( "non-linear minimization with optim(): iteration limit exceeded" )
      } else {
         stop( "non-linear minimization: 'optim' did not converge" )
      }
   }

   ## asymptotic least squares (ALS)
   ## (all coefficients may adjust to best fit of the model)
   cholMat <- triang( mindist$par, nNetput - 1 )
      # triangular Cholesky Matrix
   cVecliHessian <- vecli( t(cholMat) %*% cholMat )
      # vector of linear independent values of constrained Hessian
   coef <- estResult$coef$liCoef + estResult$coef$liCoefCov %*%
      t( hessianDeriv ) %*% solve( hessianDeriv %*% estResult$coef$liCoefCov %*%
      t( hessianDeriv ) ) %*% ( cVecliHessian - uVecliHessian )
      # vector of li indep. constrained coefficients

   ## computation of the coefficient variance covariance matrix
   if( stErMethod == "none" ) {
      coefVcov <- NULL
   } else {
      if( verbose >= 1 ) {
         cat( "Starting simulation to obtain standard errors\n" )
      }
      result$sim <- .snqProfitImposeConvexityStEr( estResult = estResult,
         rankReduction = rankReduction, start = start,
         optimMethod = optimMethod, control = control,
         stErMethod = stErMethod, nRep = nRep, verbose = verbose )
      coefVcov <- result$sim$coefVcov
   }

   ## results of constrained model
   result$pMeans <- estResult$pMeans
   result$qMeans <- estResult$qMeans
   result$fMeans <- estResult$fMeans
   result$mindist <- mindist
   result$coef <- snqProfitCoef( coef, nNetput, nFix, form = estResult$form,
      quantNames = names( estResult$qMeans ), priceNames = names( estResult$pMeans ),
      fixNames = names( estResult$fMeans ), coefCov = coefVcov,
      df = estResult$est$df )
      # constrained coefficients
   result$fitted <- snqProfitCalc( priceNames, fixNames, data = estResult$data,
      weights = estResult$weights, scalingFactors = estResult$scalingFactors,
      coef = result$coef, form = estResult$form, quantNames = quantNames )
   result$residuals <- data.frame( nr = c( 1:nObs ) )
   for( i in 1:nNetput ) {
      result$residuals[[ quantNames[ i ] ]] <- estResult$data[[ quantNames[ i ] ]] /
         estResult$scalingFactors[ i ] - result$fitted[ , i ]
   }
   if( !( "nr" %in% quantNames ) ) {
      result$residuals[[ "nr" ]] <- NULL
   }
   result$r2 <- array( NA, c( nNetput ) )
   for( i in 1:nNetput ) {
      result$r2[ i ] <- rSquared( estResult$data[[ quantNames[ i ] ]] /
         estResult$scalingFactors[ i ], result$residuals[[ quantNames[ i ] ]] )
   }
   names( result$r2 ) <- names( estResult$qMeans )

   result$hessian <- snqProfitHessian( result$coef$beta, estResult$pMeans,
      estResult$weights ) # constrained Hessian matrix
   result$ela <- snqProfitEla( result$coef$beta, estResult$pMeans, estResult$qMeans,
      estResult$weights, coefVcov = result$coef$allCoefCov,
      df = result$est$df ) # elasticities of constrained model
   if( nFix > 0 && estResult$form == 0 ) {
      result$fixEla <- snqProfitFixEla( result$coef$delta, result$coef$gamma,
         result$qMeans, result$fMeans, estResult$weights )
   }
   result$est       <- estResult$est
   result$data      <- estResult$data
   result$weights   <- estResult$weights
   result$normPrice <- estResult$normPrice
   result$convexity <- TRUE
   result$priceNames    <- estResult$priceNames
   result$quantNames    <- estResult$quantNames
   result$fixNames    <- estResult$fixNames
   result$form      <- estResult$form
   result$base      <- estResult$base
   result$method    <- estResult$method
   result$scalingFactors <- estResult$scalingFactors

   class( result ) <- "snqProfitImposeConvexity"
   return( result )
}
