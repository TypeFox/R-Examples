# fitted "frontier" values of frontier models
fitted.frontier <- function( object, asInData = FALSE, ... ) {

   if( asInData ) {
      result <- rep( NA, length( object$validObs ) )
      if( sum( object$validObs ) != nrow( object$dataTable ) ) {
         stop( "internal error: number of rows of element 'dataTable' is not",
            " equal to number of valid observations" )
      }
      for( i in 1:nrow( object$dataTable ) ) {
         result[ object$validObs ][ i ] <- object$fitted[ 
            object$dataTable[ i, 1 ], object$dataTable[ i, 2 ] ]
      }
      names( result ) <- names( object$validObs )
   } else {
      result <- object$fitted
   }
   return( result )
}
