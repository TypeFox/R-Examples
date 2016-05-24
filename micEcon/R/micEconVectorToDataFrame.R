.micEconVectorToDataFrame <- function( data ) {

   if( !is.data.frame( data ) ) {
      if( is.vector( data ) ) {
         varNames <- names( data )
         if( is.null( varNames ) ) {
            stop( "if argument 'data' is a vector,",
               " its elements must have names" )
         } else {
            data <- as.data.frame( matrix( data, nrow = 1 ) )
            names( data ) <- varNames
         }
      } else {
         stop( "argument 'data' must be either a data.frame or a vector" )
      }
   }

   return( data )
}