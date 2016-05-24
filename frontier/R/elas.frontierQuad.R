elas.frontierQuad <- function( object, data = NULL, dataLogged = TRUE,
      yObs = FALSE, ... ) {

   if( is.null( data ) ) {
      data <- eval( object$call$data )
   }

   if( !is.logical( dataLogged ) || length( dataLogged ) != 1 ) {
      stop( "argument 'dataLogged' must be a single logical value" )
   }

   if( yObs ) {
      yName <- object$yName
   } else {
      yName <- NULL
   }

   xNames <- eval( object$call$xNames )

   nExog <- length( xNames )
   nCoef <- 1 + nExog + nExog * ( nExog + 1 ) / 2

   if( dataLogged ) {
      # translog function
      result <- translogEla( xNames = xNames, data = data,
         coef = coef( object )[ 1:nCoef ],
         coefCov = vcov( object )[ 1:nCoef, 1:nCoef ],
         dataLogged = TRUE )
   } else {
      # quadratic function
      stop( "sorry, the calculation of elasticities of a quadratic function",
         " is not implemented yet" )
   }

   return( result )
}