coef.frontier <- function( object, which = "mle", extraPar = FALSE, ... ) {

   if( length( extraPar ) != 1 || !is.logical( extraPar[1] ) ) {
      stop( "argument 'extraPar' must be a single logical value" )
   }

   if( tolower( which ) == "start" ){
      result <- object$startVal
   } else if( tolower( which ) == "ols" ) {
      result <- object$olsParam
   } else if( tolower( which ) == "grid" ) {
      result <- object$gridParam
   } else if( tolower( which ) == "mle" ) {
      result <- object$mleParam
   } else {
      stop( "argument 'which' must be either 'start', 'ols', 'grid',",
         " or 'mle'" )
   }
   
   if( extraPar && !is.null( result ) ) {
      if( tolower( which ) == "ols" ) {
         warning( "extra parameters are not available for coefficients",
            " obtained by OLS" )
      } else {
         result <- c( result, 
            sigmaSqU = unname( result[ "sigmaSq" ] * result[ "gamma"] ),
            sigmaSqV = unname( result[ "sigmaSq" ] * ( 1 - result[ "gamma"] ) ),
            sigma = unname( sqrt( result[ "sigmaSq" ] ) ),
            sigmaU = unname( sqrt( result[ "sigmaSq" ] * result[ "gamma"] ) ),
            sigmaV = unname( sqrt( result[ "sigmaSq" ] * ( 1 - result[ "gamma"] ) ) ),
            lambdaSq = unname( result[ "gamma"] / ( 1 -  result[ "gamma"] ) ),
            lambda = unname( sqrt( result[ "gamma"] / ( 1 -  result[ "gamma"] ) ) ) )
         if( object$modelType == 1 && ! object$timeEffect ) {
            sigmaSqU <- result[ "sigmaSqU"]
            if( object$truncNorm & which == "mle" ) {
               alpha <- - result[ "mu" ] / sqrt( sigmaSqU )
            } else {
               alpha <- 0
            }
            denom <- 1 - pnorm( alpha )
            result <- c( result,
               varU = unname( sigmaSqU * ( 1 + alpha * dnorm( alpha ) / denom -
                     ( dnorm( alpha ) / denom )^2 ) ) )
            result <- c( result,
               sdU = unname( sqrt( result[ "varU" ] ) ),
               gammaVar = unname( result[ "varU" ] / 
                     ( result[ "varU" ] + result[ "sigmaSqV" ] ) ) )
         }
      }
   }
   
   return( result )
}