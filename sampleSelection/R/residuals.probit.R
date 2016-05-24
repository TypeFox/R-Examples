residuals.probit <- function( object, type = "deviance", ... ) {

   fitVal <- fitted( object )
   response <- model.frame( object )[ , 1 ]
   if( type == "response" ) {
      result <- response - fitVal
   } else if( type == "deviance" ) {
      result <- ifelse( response == 1,
         sqrt( -2 * log( fitVal ) ), -sqrt( -2 * log( 1 - fitVal ) ) )
   } else if( type == "pearson" ) {
      result <- ( response - fitVal ) / sqrt( fitVal * ( 1 - fitVal ) )
   } else {
      stop( "argument 'type' must be either 'deviance', 'pearson',",
         " or 'response'" )
   }
   result <- drop( result )
   
   if( !is.null( object$weights ) && type %in% c( "pearson", "deviance" ) ) {
      result <- result * sqrt( object$weights )
   }

   names( result ) <- names( fitted( object ) )
   return( result )
}
