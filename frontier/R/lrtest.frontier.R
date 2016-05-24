## Likelihood Ratio Test
lrtest.frontier <- function( object, ... ) {

   thisCall <- match.call()

   if( ! "frontier" %in% class( object ) ){
      stop( "argument 'object' must be of class 'frontier'" )
   }

   ## list of objects in ...
   objectList <- list( ... )

   if( length( objectList ) < 1 ) {
      ## if there is only one object, test for the inefficiency term
      # obtain log likelihood values (including #DFs)
      olsLogLik <- logLik( object, which = "OLS" )
      mleLogLik <- logLik( object )
      # prepare object to be returned
      result <- matrix( NA, nrow = 2, ncol = 5 )
      colnames( result ) <- c( "#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)" )
      result <- as.data.frame( result )
      # fill in #DF and log likelihood value of OLS estimation
      result[ 1, "#Df" ] <- attributes( olsLogLik )$df
      result[ 1, "LogLik" ] <- olsLogLik[ 1 ]
      # fill in #DF and log likelihood value of MLE estimation
      result[ 2, "#Df" ] <- attributes( mleLogLik )$df
      result[ 2, "LogLik" ] <- mleLogLik[ 1 ]
      # #DF of the test
      result[ 2, "Df" ] <- result[ 2, "#Df" ] - result[ 1, "#Df" ]
      # test statistic
      result[ 2, "Chisq" ] <- 2 *
         ( result[ 2, "LogLik" ] - result[ 1, "LogLik" ] )
      # P-value
      if( result[ 2, "Df" ] == 1 ) {
         result[ 2, "Pr(>Chisq)" ] <-
            0.5 * pchisq( result[ 2, "Chisq" ], 0, lower.tail = FALSE ) +
            0.5 * pchisq( result[ 2, "Chisq" ], 1, lower.tail = FALSE )
      } else if( result[ 2, "Df" ] > 1 ) {
         result[ 2, "Pr(>Chisq)" ] <-
            0.25 * pchisq( result[ 2, "Chisq" ], result[ 2, "Df" ] - 2,
               lower.tail = FALSE ) +
            0.5 * pchisq( result[ 2, "Chisq" ], result[ 2, "Df" ] - 1,
               lower.tail = FALSE ) +
            0.25 * pchisq( result[ 2, "Chisq" ], result[ 2, "Df" ],
               lower.tail = FALSE )
      } else {
         stop( "internal error: degrees of freedom of the LR test are",
            " non-positive (", result[ 2, "Df" ], ")" )
      }

      # add appropriate class
      class( result ) <- c( "anova", class( result ) )
      # add heading for this test
      attributes( result )$heading <- c( "Likelihood ratio test\n",
         paste( "Model 1: OLS (no inefficiency)\nModel 2:",
            ifelse( object$modelType == 1,
               "Error Components Frontier (ECF)",
               "Efficiency Effects Frontier (EEF)" ) ) )
   } else {
      ## if there is more than one object, compare the models

      ## test if all objects are of class "frontier"
      for( i in 1:length( objectList ) ) {
         if( ! "frontier" %in% class( objectList[[ i ]] ) ){
            stop( "all further arguments ('...') must be of class 'frontier'" )
         }
      }

      ## get and save the names of the models (objects)
      object$lrtest.frontier.name <- deparse( substitute( object ) )
      dotsNames <- as.list( thisCall )[ -1 ]
      dotsNames$object <- NULL
      for( i in 1:length( objectList ) ){
         objectList[[ i ]]$lrtest.frontier.name <- deparse( dotsNames[[ i ]] )
      }
      ## function to extract the names of the objects
      extractName <- function( object ){
         return( object$lrtest.frontier.name )
      }

      ## test if all models are of the same model type
      for( i in 1:length( objectList ) ) {
         if( object$modelType != objectList[[ i ]]$modelType ){
            stop( "all models must be of the same type",
               " but model '", extractName( object ), "' is an",
               ifelse( object$modelType == 1,
                  " 'Error Components Frontier' (ECF)",
                  " 'Efficiency Effects Frontier' (EEF)" ),
               ", while model '", extractName( objectList[[ i ]] ), "' is an",
               ifelse( objectList[[ i ]]$modelType == 1,
                  " 'Error Components Frontier' (ECF)",
                  " 'Efficiency Effects Frontier' (EEF)" ) )
         }
      }

      ## do the LR tests
      result <- do.call( lrtest.default,
         c( list( object = object ), objectList, list( name = extractName ) ) )

      for( i in 2:nrow( result ) ){
         if( ( result[ i, "#Df" ] - result[ i - 1, "#Df" ] ) *
               ( result[ i, "LogLik" ] - result[ i - 1, "LogLik" ] ) < 0 ) {
            if( result[ i, "LogLik" ] > result[ i - 1, "LogLik" ] ) {
               compareLikelihood <- "larger"
               compareDf <- "less"
            } else {
               compareLikelihood <- "smaller"
               compareDf <- "more"
            }
            warning( "model '", i, "' has a ", compareLikelihood,
               " log-likelihood value than the ", compareDf,
               " restricted model '", i - 1, "'" )
         }
      }
   }

   return( result )
}
