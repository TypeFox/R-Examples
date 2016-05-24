elas.translogEst <- function( object, data = NULL, dataLogged = NULL, ... ) {

   if( is.null( data ) ) {
      data <- eval( object$call$data )
   }

   if( is.null( dataLogged ) ) {
      dataLogged <- object$dataLogged
   }

   result <- translogEla( xNames = object$xNames,
      data = data, coef = coef( object ), coefCov = vcov( object ),
      dataLogged = dataLogged )

   return( result )
}