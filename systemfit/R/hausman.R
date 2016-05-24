## this function returns test statistic for
## the hausman test
# The m-statistic is then distributed with k degrees of freedom, where k
# is the rank of the matrix .A generalized inverse is used, as
# recommended by Hausman (1982).

hausman.systemfit <- function( results2sls, results3sls ) {

   result <- list()

   if( is.null( results2sls$restrict.regMat ) ) {
      result$q <- coef( results2sls ) - coef( results3sls )
      result$qVar <- vcov( results2sls ) - vcov( results3sls )
   } else {
      result$q <- coef( results2sls, modified.regMat = TRUE ) -
         coef( results3sls, modified.regMat = TRUE )
      result$qVar <- vcov( results2sls, modified.regMat = TRUE ) -
         vcov( results3sls, modified.regMat = TRUE )
   }

#    if( min( eigen( hausman$qVar )$values ) < 0 ) {
#       warning( "the matrix V is not 'positive definite'" )
#    }

   result$statistic <- crossprod( result$q, solve( result$qVar, result$q ) )
   names( result$statistic ) <- "Hausman"
   result$parameter <- nrow( result$qVar )
   names( result$parameter ) <- "df"
   result$p.value <- 1 - pchisq( result$statistic, result$parameter )
   result$method = paste( "Hausman specification test for consistency of",
      "the 3SLS estimation" )
   if( "data" %in% names( results2sls$call ) ) {
      result$data.name <- results2sls$call$data
   } else {
      result$data.name <- "unknown"
   }
   class( result ) <- "htest"
   return( result )
}
