## ===== calculation of fixed factor elasticities ===
snqProfitFixEla <- function( delta, gamma, quant, fix, weights,
      scalingFactors = rep( 1, length( weights ) ) ) {

   if( ncol( delta ) > 1 ) {
      if( nrow( gamma ) != ncol( gamma ) ) {
         stop( "argument 'gamma' must be a quadratic matrix" )
      }
   }
   if( ( ncol( delta ) == 1 && length( gamma ) != 1 ) ||
         ( ncol( delta ) > 1 && ncol( gamma ) != ncol( delta ) ) ) {
      stop( "the number of columns of argument 'delta' must be equal to",
         " the number of rows and the number of columns of argument 'gamma'" )
   }
   if( length( quant ) != nrow( delta ) ) {
      stop( "the length of arguments 'quant' must be equal to the number of",
         " columns of argument 'quant'" )
   }
   if( length( fix ) != ncol( delta ) ) {
      stop( "the length of arguments 'fix' must be equal to the number of",
         " columns of argument 'delta'" )
   }
   if( length( quant ) != length( weights ) ) {
      stop( "arguments 'quant' and 'weights' must have the same length" )
   }
   nNetput  <- length( quant )
   nFix     <- length( fix )
   quant    <- unlist( quant ) / scalingFactors
   fix      <- unlist( fix )

   ela <- ( delta + weights %*% t( fix ) %*% gamma ) *
      array( 1, nNetput ) %*% t( fix ) /
      quant %*% t( array( 1, nFix ) )

   return( ela )
}
