elas.quadFuncEst <- function( object, data = NULL, yObs = FALSE, ... ) {

   if( is.null( data ) ) {
      data <- eval( object$call$data )
   }

   if( yObs ) {
      yName = object$yName
   } else {
      yName <- NULL
   }

   result <- quadFuncEla( xNames = object$xNames,
      data = data, coef = coef( object ), yName = yName,
      shifterNames = object$shifterNames,
      homWeights = object$homWeights )

   return( result )
}