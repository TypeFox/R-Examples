snqProfitEst <- function( priceNames, quantNames, fixNames = NULL,
   instNames = NULL, data,  form = 0, base = 1, scalingFactors = NULL,
   weights = snqProfitWeights( priceNames, quantNames, data, "DW92", base = base ),
   method = ifelse( is.null( instNames ), "SUR", "3SLS" ), ... ) {

   checkNames( c( priceNames, quantNames, fixNames, instNames ), names( data ) )

   if( length( quantNames ) != length( priceNames ) ) {
      stop( "arguments 'quantNames' and 'priceNames' must have the same length" )
   }
   if( length( priceNames ) < 2 ) {
      stop( "you must specify at least 2 netputs" )
   }
   if( length( priceNames ) != length( weights ) ) {
      stop( "arguments 'priceNames' and 'weights' must have the same length" )
   }
   if( min( weights ) < 0 ) {
      warning( "At least one weight of the prices for normalization",
         " (argument 'weights') is negative. Thus, in this case positive",
         " semidefiniteness of the 'beta' matrix does not ensure",
         " a convex profit function." )
   }
   if( !is.null( scalingFactors ) ) {
      if( length( scalingFactors ) != length( priceNames ) ) {
         stop( "arguments 'priceNames' and 'scalingFactors' must have the",
            " same length" )
      }
      if( base != 1 ) {
         warning( "argument 'base' is ignored because argument",
            " 'scalingFactors' is provided" )
      }
   }

   nNetput <- length( quantNames )  # number of netputs
   nFix    <- length( fixNames )  # number of fixed inputs
   nIV     <- length( instNames )  # number of fixed inputs
   nObs    <- nrow( data )      # number of observations

   if( form == 0 ) {
      nCoef   <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         ( nFix + 1 ) * nFix/2  #number of coefficients
   } else if( form == 1 ) {
      nCoef   <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         nNetput * ( nFix + 1 ) * nFix/2  #number of coefficients
   } else {
      stop( "argument 'form' must be either 0 or 1" )
   }

   result  <- list()

   ## scaling factors
   if( is.null( scalingFactors ) ) {
      scalingFactors <- rep( 1, nNetput )
      if( !is.null( base ) ) {
         for( i in 1:nNetput ) {
            scalingFactors[ i ] <- 1 / mean( data[[ priceNames[ i ] ]][ base ] )
         }
      }
   }

   ## mean Values
   result$pMeans <- array( NA, nNetput )
   result$qMeans <- array( NA, nNetput )
   for( i in 1:nNetput ) {
      result$pMeans[ i ] <- mean( data[[ priceNames[ i ] ]] ) * scalingFactors[ i ]
      result$qMeans[ i ] <- mean( data[[ quantNames[ i ] ]] ) / scalingFactors[ i ]
   }
   names( result$pMeans ) <- priceNames
   names( result$qMeans ) <- quantNames
   if( nFix > 0 ) {
      result$fMeans <- array( NA, nFix )
      for( i in 1:nFix ) {
         result$fMeans[ i ] <- mean( data[[ fixNames[ i ] ]] )
      }
      names( result$fMeans ) <- fixNames
   }

   ## instrumental variables
   if( nIV == 0 ) {
      inst <- NULL
   } else {
      inst <- as.formula( paste( "~", paste( paste( "iv", c( 1:nIV ), sep="" ),
         collapse = " + " ) ) )
   }

   ## prepare and estimate the model
   modelData <- .snqProfitModelData( data = data, weights = weights,
      priceNames = priceNames, quantNames = quantNames, fixNames = fixNames, instNames = instNames,
      form = form, netputScale = scalingFactors, fixedScale = result$fMeans )
   system <- snqProfitSystem( nNetput, nFix )    # equation system
   restrict <- snqProfitRestrict( nNetput, nFix, form )    # restrictions
   result$est <- systemfit( formula = system, method = method, data = modelData,
      restrict.regMat = restrict, inst = inst, ... )
   result$coef <- snqProfitCoef( coef = coef( result$est, modified.regMat = TRUE ),
      nNetput = nNetput, nFix = nFix, form = form,
      coefCov = vcov( result$est, modified.regMat = TRUE ),
      df = nNetput * nObs - nCoef,
      quantNames = quantNames, priceNames = priceNames, fixNames = fixNames )
      # estimated coefficients
   result$coef <- .snqProfitRescaleCoef( result$coef, nNetput, result$fMeans,
      form )

   result$fitted <- data.frame( profit0 = rep( 0, nObs ) )
   result$residuals <- data.frame( profit0 = rep( 0, nObs ) )
   for( i in 1:nNetput ) {
      result$fitted[[ quantNames[ i ] ]] <- fitted( result$est$eq[[ i ]] )
      result$fitted[[ "profit0" ]] <- result$fitted[[ "profit0" ]] +
         result$fitted[[ quantNames[ i ] ]] * data[[ priceNames[ i ] ]] *
         scalingFactors[ i ]
      result$residuals[[ quantNames[ i ] ]] <- data[[ quantNames[ i ] ]] /
         scalingFactors[ i ] - result$fitted[[ quantNames[ i ] ]]
   }
   result$fitted[[ "profit" ]] <- result$fitted[[ "profit0" ]]
   result$fitted[[ "profit0" ]] <- NULL
   result$residuals[[ "profit" ]] <- modelData[[ "profit" ]] -
      result$fitted[[ "profit" ]]
   result$residuals[[ "profit0" ]] <- NULL

   result$r2 <- array( NA, c( nNetput + 1 ) )
   for( i in 1:nNetput ) {
      # result$r2[ i ] <- summary( result$est$eq[[ i ]] )$r.squared
      result$r2[ i ] <- rSquared( data[[ quantNames[ i ] ]] / scalingFactors[ i ],
         result$residuals[[ quantNames[ i ] ]] )
   }
   result$r2[ nNetput + 1 ] <- rSquared( modelData[[ "profit" ]],
      result$residuals[[ "profit" ]] )
   names( result$r2 ) <- c( quantNames, "profit" )

   result$hessian <- snqProfitHessian( result$coef$beta, result$pMeans, weights )
      # Hessian matrix
   result$ela <- snqProfitEla( result$coef$beta, result$pMeans,
      result$qMeans, weights, coefVcov = result$coef$allCoefCov,
      df = df.residual( result$est ) )   # estimated elasticities
   if( nFix > 0 && form == 0 ) {
      result$fixEla <- snqProfitFixEla( result$coef$delta, result$coef$gamma,
         result$qMeans, result$fMeans, weights )
   }

   result$data     <- data
   result$weights  <- weights
   names( result$weights ) <- priceNames
   result$normPrice <- modelData$normPrice
   if( nNetput > 2 ){
      result$convexity  <- semidefiniteness( result$hessian[
         1:( nNetput - 1 ), 1:( nNetput - 1 ) ], positive = TRUE )
   } else if( nNetput == 2 ){
      result$convexity  <- result$hessian[ 1, 1 ] >= 0
   }
   result$priceNames  <- priceNames
   result$quantNames  <- quantNames
   result$fixNames  <- fixNames
   result$instNames <- instNames
   result$form    <- form
   result$base    <- base
   result$method  <- method
   result$scalingFactors <- scalingFactors
   names( result$scalingFactors ) <- priceNames

   class( result )  <- "snqProfitEst"
   return( result )
}
