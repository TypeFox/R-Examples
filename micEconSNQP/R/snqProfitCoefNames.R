snqProfitCoefNames <- function( nNetput, nFix, form = 0, all=FALSE ) {
   names <- NULL
   for( i in 1:nNetput ) {
      names <- c( names, paste( "alpha" , i ) )
   }
   for( i in 1:nNetput ) {
      for( j in 1:nNetput ) {
         if( all || ( j >= i && j < nNetput && i <nNetput ) ) {
            names <- c( names, paste( "beta", i, j ) )
         }
      }
   }
   if( nFix > 0 ) {
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            names <- c( names, paste( "delta", i, j ) )
         }
      }
      if( form == 0 ) {
         for( i in 1:nFix ) {
            for( j in 1:nFix ) {
               if( all || j>= i ) {
                  names <- c( names, paste( "gamma", i, j ) )
               }
            }
         }
      } else if( form == 1 ) {
         for( n in 1:nNetput ) {
            for( i in 1:nFix ) {
               for( j in 1:nFix ) {
                  if( all || j>= i ) {
                     names <- c( names, paste( "gamma", n, i, j ) )
                  }
               }
            }
         }
      } else {
         stop( "argument 'form' must be either 0 or 1" )
      }
   }
   return( names )
}
