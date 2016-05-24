.snqProfitImposeConvexityStEr <- function( estResult, rankReduction,
   start, optimMethod, control, stErMethod, nRep, verbose ) {

   ## computation of the coefficient variance covariance matrix
   nObs <- nrow( estResult$data )
   nCoef <- length( estResult$coef$liCoef )
   nAllCoef <- length( estResult$coef$allCoef )
   nNetput  <- length( estResult$pMeans )
   nFix     <- length( estResult$fMeans )
   sim <- list()
   if( stErMethod == "jackknife" ) {
      nRep <- nObs
   }
   sim$coef    <- matrix( NA, nCoef, nRep )
   sim$allCoef <- matrix( NA, nAllCoef, nRep )
   sim$status  <- rep( NA, nRep )
   sim$results <- list()
   if( stErMethod %in% c( "jackknife", "resample" ) ) {
      for( repNo in 1:nRep ) {
         if( verbose >= 2 ) {
            cat( repNo, " / ", nRep, sep = "" )
         }
         if( stErMethod == "jackknife" ) {
            simData <- estResult$data[ -repNo, ]
         } else if( stErMethod == "resample" ) {
            simData <- estResult$data[ ceiling( runif( nObs, 0, nObs ) ), ]
         }
         simResult <- try( snqProfitEst( priceNames = estResult$priceNames,
            quantNames = estResult$quantNames, fixNames = estResult$fixNames,
            instNames = estResult$instNames, data = simData, form = estResult$form,
            scalingFactors = estResult$scalingFactors,
            weights = estResult$weights, method = estResult$method ),
            silent = TRUE )
         if( class( simResult) == "try-error" ) {
            sim$status[ repNo ] <- 888
         } else {
            if( simResult$convexity ) {
               sim$status[ repNo ] <- -1
            } else {
               simResult <- try( snqProfitImposeConvexity( simResult,
                  rankReduction = rankReduction, start = start,
                  optimMethod = optimMethod, control = control ),
                  silent = TRUE )
               if( class( simResult) == "try-error" ) {
                  sim$status[ repNo ] <- 999
               } else {
                  sim$status[ repNo ] <- simResult$mindist$convergence
               }
            }
         }
         if( sim$status[ repNo ] <= 0 ) {
            sim$coef[ , repNo ] <- simResult$coef$liCoef
            sim$allCoef[ , repNo ] <- simResult$coef$allCoef
            # sim$results[[ repNo ]] <- simResult
         }
         if( verbose >= 3 ) {
            cat( ", status: ", sim$status[ repNo ], "\n", sep = "" )
         }
      }
   } else if( stErMethod == "coefSim" ) {
      fakeResult <- list()
      class( fakeResult )  <- "snqProfitEst"
      fakeResult$qMeans    <- estResult$qMeans
      fakeResult$pMeans    <- estResult$pMeans
      fakeResult$fMeans    <- estResult$fMeans
      fakeResult$weights   <- estResult$weights
      fakeResult$form      <- estResult$form
      fakeResult$data      <- estResult$data
      fakeResult$normPrice <- estResult$normPrice
      fakeResult$scalingFactors <- estResult$scalingFactors
      for( repNo in 1:nRep ) {
         if( verbose >= 2 ) {
            cat( repNo, " / ", nRep, sep = "" )
         }
         liCoef <- mvrnorm( mu = estResult$coef$liCoef,
            Sigma = estResult$coef$liCoefCov )
         fakeResult$coef <- snqProfitCoef( liCoef, nNetput, nFix,
            form = estResult$form, coefCov = estResult$coef$liCoefCov )
         fakeResult$hessian <- snqProfitHessian( fakeResult$coef$beta,
            fakeResult$pMeans, fakeResult$weights )
         fakeResult$convexity <- semidefiniteness( fakeResult$hessian[
            1:( nNetput - 1 ), 1:( nNetput - 1 ) ], positive = TRUE )
         if( fakeResult$convexity ) {
            simResult <- fakeResult
            sim$status[ repNo ] <- -1
         } else {
            simResult <- try( snqProfitImposeConvexity( fakeResult,
               rankReduction = rankReduction, start = start,
               optimMethod = optimMethod, control = control ),
               silent = ( verbose < 2 ) )
            if( class( simResult) == "try-error" ) {
               sim$status[ repNo ] <- 999
            } else {
               sim$status[ repNo ] <- simResult$mindist$convergence
            }
         }
         if( sim$status[ repNo ] <= 0 ) {
            sim$coef[ , repNo ]    <- simResult$coef$liCoef
            sim$allCoef[ , repNo ] <- simResult$coef$allCoef
            # sim$results[[ repNo ]] <- simResult
         }
         if( verbose >= 2 ) {
            cat( ", status: ", sim$status[ repNo ], "\n", sep = "" )
         }
      }
   }
   nValidRep       <- sum( sim$status <= 0 )
   sim$coef        <- sim$coef[ , sim$status <= 0 ]
   sim$allCoef     <- sim$allCoef[ , sim$status <= 0 ]
   sim$coefMean    <- rowMeans( sim$coef )
   sim$allCoefMean <- rowMeans( sim$allCoef )
   sim$coefDev     <- sim$coef - matrix( sim$coefMean, nCoef, nValidRep )
   sim$allCoefDev  <- sim$allCoef - matrix( sim$allCoefMean, nAllCoef,
      nValidRep )
   sim$coefVcov    <- sim$coefDev    %*% t( sim$coefDev )
   sim$allCoefVcov <- sim$allCoefDev %*% t( sim$allCoefDev )
   if( stErMethod == "jackknife" ) {
      sim$coefVcov    <- ( ( nValidRep - 1 ) / nValidRep ) * sim$coefVcov
      sim$allCoefVcov <- ( ( nValidRep - 1 ) / nValidRep ) *
         sim$allCoefVcov
   } else if( stErMethod %in% c( "resample", "coefSim" ) ) {
      sim$coefVcov    <- sim$coefVcov    / nValidRep
      sim$allCoefVcov <- sim$allCoefVcov / nValidRep
   }
   return( sim )
}
