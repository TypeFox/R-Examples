print.summary.mvProbitMargEff <- function( x, digits = 4, ... ) {

   if( length( digits ) == 5 ) {
     for( i in 1:4 ) {
       x[ , i ] <- round( x[ , i ], digits = digits[ i + 1 ] )
     }
     digits <- digits[1]
   }
  
   printCoefmat( x, digits = digits )

   invisible( x )
}

