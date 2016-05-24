cesEstGridRho <- function( rho1Values, rho2Values, rhoValues, returnAll, ... )  {

   # some tests
   if( is.null( rho1Values ) && is.null( rho2Values ) && is.null( rhoValues ) ) {
      stop( "either argument 'rho1Values', 'rho2Values', or 'rhoValues'",
         " must be non-NULL" )
   }
   if( !is.null( rho1Values ) ) {
      if( !is.numeric( rho1Values ) ) {
         stop( "the rho_1s specified in argument 'rho1Values'",
            " must be numeric" )
      } else if(  min( rho1Values ) < -1 ) {
         stop( "the rho_1s specified in argument 'rho1Values'",
            " must not be smaller than '-1'" )
      }
   }
   if( !is.null( rho2Values ) ) {
      if( !is.numeric( rho2Values ) ) {
         stop( "the rho_2s specified in argument 'rho2Values'",
            " must be numeric" )
      } else if(  min( rho2Values ) < -1 ) {
         stop( "the rho_2s specified in argument 'rho2Values'",
            " must not be smaller than '-1'" )
      }
   }
   if( !is.null( rhoValues ) ) {
      if( !is.numeric( rhoValues ) ) {
         stop( "the rhos specified in argument 'rhoValues'",
            " must be numeric" )
      } else if( min( rhoValues ) < -1 ) {
         stop( "the rhos specified in argument 'rhoValues'",
            " must not be smaller than '-1'" )
      }
   }

   # list that should contain each estimation result
   allResults <- list()

   # summary results for each estimation (with different fixed rhos)
   if( !is.null( rho1Values ) && is.null( rho2Values ) && is.null( rhoValues ) ) {
      sumResults <- data.frame( rho1 = rho1Values )
   } else if( is.null( rho1Values ) && !is.null( rho2Values ) && is.null( rhoValues ) ) {
      sumResults <- data.frame( rho2 = rho2Values )
   } else if( is.null( rho1Values ) && is.null( rho2Values ) && !is.null( rhoValues ) ) {
      sumResults <- data.frame( rho = rhoValues )
   } else if( !is.null( rho1Values ) && !is.null( rho2Values ) && is.null( rhoValues ) ) {
      sumResults <- expand.grid( rho1 = rho1Values, rho2 = rho2Values )
   } else if( !is.null( rho1Values ) && is.null( rho2Values ) && !is.null( rhoValues ) ) {
      sumResults <- expand.grid( rho1 = rho1Values, rho = rhoValues )
   } else if( is.null( rho1Values ) && !is.null( rho2Values ) && !is.null( rhoValues ) ) {
      sumResults <- expand.grid( rho2 = rho2Values, rho = rhoValues )
   } else if( !is.null( rho1Values ) && !is.null( rho2Values ) && !is.null( rhoValues ) ) {
      sumResults <- expand.grid( rho1 = rho1Values, rho2 = rho2Values, rho = rhoValues )
   }
   sumResults$rss <- NA
   sumResults$convergence <- NA

   # estimate the CES for each pre-defined rho
   for( i in 1:nrow( sumResults ) ) {
      allResults[[ i ]] <- cesEst( rho1 = sumResults[[ "rho1" ]][ i ], 
         rho2 = sumResults[[ "rho2" ]][ i ], 
         rho = sumResults[[ "rho" ]][ i ], ... )
      sumResults$rss[ i ] <- allResults[[ i ]]$rss
      if( !is.null( allResults[[ i ]]$convergence ) ) {
         sumResults$convergence[ i ] <- allResults[[ i ]]$convergence
      }
   }

   # returned object: the estimation results with the lowest RSS
   result <- allResults[[ which.min( sumResults$rss ) ]]

   # add the summary results of each estimation
   result$allRhoSum <- sumResults

   # add full results of each estimation
   if( returnAll ) {
      result$allRhoFull <- allResults
   }
   
   if( !is.null( rho1Values ) && !is.null( rho2Values ) && is.null( rhoValues ) ) {
      result$rssArray <- matrix( result$allRhoSum$rss, 
         nrow = length( rho1Values ), ncol = length( rho2Values ),
         byrow  = FALSE )
      rownames( result$rssArray ) <- rho1Values
      colnames( result$rssArray ) <- rho2Values
   } else if( !is.null( rho1Values ) && is.null( rho2Values ) && !is.null( rhoValues ) ) {
      result$rssArray <- matrix( result$allRhoSum$rss, 
         nrow = length( rho1Values ), ncol = length( rhoValues ),
         byrow  = FALSE )
      rownames( result$rssArray ) <- rho1Values
      colnames( result$rssArray ) <- rhoValues
   } else if( is.null( rho1Values ) && !is.null( rho2Values ) && !is.null( rhoValues ) ) {
      result$rssArray <- matrix( result$allRhoSum$rss, 
         nrow = length( rho2Values ), ncol = length( rhoValues ),
         byrow  = FALSE )
      rownames( result$rssArray ) <- rho2Values
      colnames( result$rssArray ) <- rhoValues
   } else if( !is.null( rho1Values ) && !is.null( rho2Values ) && !is.null( rhoValues ) ) {
      result$rssArray <- array( result$allRhoSum$rss, 
         dim = c( length( rho1Values ), length( rho2Values ), length( rhoValues ) ),
         dimnames <- list( rho1Values, rho2Values, rhoValues ) )
   }


   # add values used for rho_1
   result$rho1Values <- rho1Values

   # add values used for rho_2
   result$rho2Values <- rho2Values

   # add values used for rho
   result$rhoValues <- rhoValues

   result$call <- match.call()
   return( result )
}
