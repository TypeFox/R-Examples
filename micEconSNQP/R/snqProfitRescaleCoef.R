.snqProfitRescaleCoef <- function( myCoef, nNetput, fixedQuant, form ) {

   nFix <- length( fixedQuant )
   if( nFix == 0 ) {
      return( myCoef )
   }
   
   # delta[ i, j ]
   for( i in 1:nNetput ) {
      for( j in 1:nFix ) {
         myCoef$delta[ i, j ] <- myCoef$delta[ i, j ] /
            fixedQuant[ j ]
         # all myCoefficients
         k <- nNetput + nNetput^2 + ( i - 1 ) * nFix + j
         myCoef$allCoef[ k ] <- myCoef$allCoef[ k ] /
            fixedQuant[ j ]
         myCoef$allCoefCov[ k, ] <- myCoef$allCoefCov[ k, ] /
            fixedQuant[ j ]
         myCoef$allCoefCov[ , k ] <- myCoef$allCoefCov[ , k ] /
            fixedQuant[ j ]
         myCoef$stats[ k, c( 1, 2 ) ] <- myCoef$stats[ k, c( 1, 2 ) ] /
            fixedQuant[ j ]
         # linear independent myCoefficients
         k <- nNetput + ( nNetput * ( nNetput - 1 ) ) / 2 +
            ( i - 1 ) * nFix + j
         myCoef$liCoef[ k ] <- myCoef$liCoef[ k ] /
            fixedQuant[ j ]
         myCoef$liCoefCov[ k, ] <- myCoef$liCoefCov[ k, ] /
            fixedQuant[ j ]
         myCoef$liCoefCov[ , k ] <- myCoef$liCoefCov[ , k ] /
            fixedQuant[ j ]
      }
   }
   if( form == 0 ) {
      # gamma[ i, j ]
      for( i in 1:nFix ) {
         for( j in 1:nFix ) {
            # gamma[ i, j ]
            myCoef$gamma[ i, j ] <- myCoef$gamma[ i, j ] /
               ( fixedQuant[ i ] * fixedQuant[ j ] )
            # all myCoefficients
            k <- nNetput + nNetput^2 + nNetput * nFix + ( i - 1 ) * nFix + j
            myCoef$allCoef[ k ] <- myCoef$allCoef[ k ] /
               ( fixedQuant[ i ] * fixedQuant[ j ] )
            myCoef$allCoefCov[ k, ] <- myCoef$allCoefCov[ k, ] /
               ( fixedQuant[ i ] * fixedQuant[ j ] )
            myCoef$allCoefCov[ , k ] <- myCoef$allCoefCov[ , k ] /
               ( fixedQuant[ i ] * fixedQuant[ j ] )
            myCoef$stats[ k, c( 1, 2 ) ] <- myCoef$stats[ k, c( 1, 2 ) ] /
               ( fixedQuant[ i ] * fixedQuant[ j ] )
            # linear independent myCoefficients
            if( j >= i ) {
               k <- nNetput + nNetput * (nNetput - 1 ) / 2 +
                  nNetput * nFix + veclipos( i, j, nFix )
               myCoef$liCoef[ k ] <- myCoef$liCoef[ k ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
               myCoef$liCoefCov[ k, ] <- myCoef$liCoefCov[ k, ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
               myCoef$liCoefCov[ , k ] <- myCoef$liCoefCov[ , k ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
            }
         }
      }
   } else {
      # gamma[ n, i, j ]
      for( n in 1:nNetput ) {
         for( i in 1:nFix ) {
            for( j in 1:nFix ) {
               # delta[ i, j ]
               myCoef$gamma[ n, i, j ] <- myCoef$gamma[ n, i, j ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
               # all myCoefficients
               k <- nNetput + nNetput^2 + nNetput * nFix +
                  ( n - 1 ) * nFix^2 + ( i - 1 ) * nFix + j
               myCoef$allCoef[ k ] <- myCoef$allCoef[ k ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
               myCoef$allCoefCov[ k, ] <- myCoef$allCoefCov[ k, ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
               myCoef$allCoefCov[ , k ] <- myCoef$allCoefCov[ , k ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
               myCoef$stats[ k, c( 1, 2 ) ] <- myCoef$stats[ k, c( 1, 2 ) ] /
                  ( fixedQuant[ i ] * fixedQuant[ j ] )
               # linear independent myCoefficients
               if( j >= i ) {
                  k <- nNetput + ( nNetput * (nNetput - 1 ) ) / 2 +
                     nNetput * nFix + ( n - 1 ) * nFix * ( nFix + 1 ) / 2 +
                     veclipos( i, j, nFix )
                  myCoef$liCoef[ k ] <- myCoef$liCoef[ k ] /
                     ( fixedQuant[ i ] * fixedQuant[ j ] )
                  myCoef$liCoefCov[ k, ] <- myCoef$liCoefCov[ k, ] /
                     ( fixedQuant[ i ] * fixedQuant[ j ] )
                  myCoef$liCoefCov[ , k ] <- myCoef$liCoefCov[ , k ] /
                     ( fixedQuant[ i ] * fixedQuant[ j ] )
               }
            }
         }
      }
   }
   return( myCoef )
}
