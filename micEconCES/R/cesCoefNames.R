cesCoefNames <- function( nExog, vrs, returnRho1 = TRUE, returnRho2 = TRUE, 
      returnRho = TRUE, nested = FALSE, withTime = FALSE ) {

   result <- "gamma"
   if( withTime ) {
      result <- c( result, "lambda" )
   }
   if( nExog == 2 ) {
      result <- c( result, "delta" )
   } else if( !nested ) {
      result <- c( result, paste( "delta", 1:nExog, sep = "_" ) )
   } else if( nested && nExog == 3 ) {
      result <- c( result, "delta_1", "delta" )
   } else if( nested && nExog == 4 ) {
      result <- c( result, "delta_1", "delta_2", "delta" )
   } else {
      stop( "internal error: non-supported arguments to cesCoefNames()" )
   }
   if( returnRho1 && nested ) {
      result <- c( result, "rho_1" )
   }
   if( returnRho2 && nested && nExog == 4 ) {
      result <- c( result, "rho_2" )
   }
   if( returnRho ) {
      result <- c( result, "rho" )
   }
   if( vrs ) {
      result <- c( result, "nu" )
   }
   return( result )
}
