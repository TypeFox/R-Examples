summary.mvProbitMargEff <- function( object, ... ) {

   nObs <- nrow( object )
   nMargEff <- ncol( object )

   result <- matrix( NA, nrow = nObs * nMargEff, ncol = 4 )
   if( nObs == 1 ) {
      rownames( result ) <- colnames( object )
   } else {
      rownames( result ) <- paste(
         rep( rownames( object ), each = nMargEff ), ": ", 
         rep( colnames( object ), nObs ), sep = "" )
   }
   colnames( result ) <- c( "Estimate", "Std. error", "z value", "Pr(> z)" )

   result[ , 1 ] <- c( t( as.matrix( object ) ) )

   margEffVCov <- attr( object, "vcov" )
   if( !is.null( margEffVCov ) ) {
      for( i in 1:nObs ) {
         result[ ( (i-1) * nMargEff + 1 ):( i * nMargEff ), 2 ] <-
            sqrt( diag( margEffVCov[ i, , ] ) )
      }
      result[ , 3 ] <- result[ , 1 ] / result[ , 2 ]
      result[ , 4 ] <- 2 * pnorm( -abs( result[ , 3 ] ) )
   }

   class( result ) <- c( "summary.mvProbitMargEff", class( result ) )

   return( result )
}

