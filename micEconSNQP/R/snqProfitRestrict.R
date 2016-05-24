## ---- snqProfit: restriction matrix --------
snqProfitRestrict <- function( nNetput, nFix, form = 0 ) {

   nxe <- 1 + nNetput + nNetput * nNetput + nFix + nFix * nFix
      #number of exogenous variables per equation

   if( form == 0 ) {
      nCoef <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         ( nFix + 1 ) * nFix/2  #number of coefficients
   } else if( form == 1 ) {
      nCoef <- nNetput + nNetput * ( nNetput - 1 ) / 2 + nNetput * nFix +
         nNetput * ( nFix + 1 ) * nFix/2  #number of coefficients
   } else {
      stop( "argument 'form' must be either 0 or 1" )
   }

   MR <- array( 0, c( nNetput * nxe, nCoef ) ) # restriction matrix

   for( i in 1:nNetput) {
      MR[ 1 + ( i - 1 ) * nxe, i ] <- 1       #alphas
   }
   for(i in 1:( nNetput - 1 ) ) {
      for( j in 1:( nNetput - 1 ) ) {
         MR[ 1 + (i-1) * nxe + j, nNetput + veclipos( i, j, nNetput - 1 ) ] <- 1
            #betas
         MR[ 1 + (i-1) * nxe + nNetput, nNetput +
            veclipos( i, j, nNetput - 1 ) ] <- -1   #beta( ,nNetput)
      }
   }
   for( j in 1:( nNetput - 1 ) ) {
      for( k in 1:( nNetput - 1 ) ) {
         MR[ 1 + ( nNetput - 1 ) * nxe + j, nNetput +
            veclipos( j, k, nNetput - 1 ) ] <- -1    #beta(nNetput, )
      }
   }
   for( j in 1:( nNetput - 1 ) ) {
      for( k in 1:( nNetput - 1 ) ) {            #beta(nNetput,nNetput)
         MR[ 1 + ( nNetput - 1 ) * nxe + nNetput,
            nNetput + veclipos( j, k, nNetput - 1 ) ] <-
            MR[ 1 +( nNetput - 1 ) * nxe + nNetput,
            nNetput + veclipos( j, k, nNetput - 1 ) ] + 1
      }
   }
   for( i in 1:nNetput ) {
      for( j in 1:( nNetput - 1 ) ) {
         for( k in 1:( nNetput - 1 ) ) {
            MR[ (i-1) * nxe + 1 + nNetput + (j-1) * nNetput + k,
               nNetput + veclipos( j, k, nNetput - 1 ) ] <- 1   #betas  (2.Term)
         }
         for( k in 1:( nNetput - 1 ) ) {
            MR[ (i-1) * nxe + 1 + nNetput + (j-1) * nNetput + nNetput,
               nNetput + veclipos( j, k, nNetput - 1 ) ] <- -1
               #beta( ,nNetput)  (2.Term)
         }
      }
      for( j in 1:( nNetput - 1 ) ) {
         for( k in 1:( nNetput - 1 ) ) {
            MR[ (i-1) * nxe + 1 + nNetput + ( nNetput - 1 ) * nNetput + j,
               nNetput + veclipos( j, k, nNetput - 1 ) ] <- -1    #beta(nNetput, )
         }
      }
      for( j in 1:( nNetput - 1 ) ) {
         for( k in 1:( nNetput - 1 ) ) {                   #beta(nNetput,nNetput)
            MR[ (i-1) * nxe + 1 + nNetput + nNetput * nNetput,
               nNetput + veclipos( j, k, nNetput - 1 ) ] <-
               MR[ (i-1) * nxe + 1 + nNetput + nNetput * nNetput,
               nNetput + veclipos( j, k, nNetput - 1 ) ] + 1
         }
      }
   }
   if( nFix > 0 ) {
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            MR[ (i-1) * nxe + 1 + nNetput + nNetput * nNetput + j,
               nNetput + nNetput * ( nNetput - 1 ) / 2 + (i-1) * nFix + j ] <- 1
               #deltas  (3.Term)
         }
      }
      if( form == 0 ) {
         for( i in 1:nNetput ) {
            for( j in 1:nFix ) {
               for( k in 1:nFix ) {                         #gammas  (4.Term)
                  MR[ (i-1) * nxe + 1 + nNetput + nNetput * nNetput + nFix +
                     (j-1) * nFix + k, nNetput + nNetput * ( nNetput - 1 ) / 2 +
                     nNetput * nFix + veclipos( j, k, nFix ) ] <- 1
               }
            }
         }
      } else if( form == 1 ) {
         for( i in 1:nNetput ) {
            for( j in 1:nFix ) {
               for( k in 1:nFix ) {                         #gammas  (4.Term)
                  MR[ (i-1) * nxe + 1 + nNetput + nNetput * nNetput + nFix +
                     (j-1) * nFix + k, nNetput + nNetput * ( nNetput - 1 ) / 2 +
                     nNetput * nFix + (i-1) * ( nFix + 1 ) * nFix/2 +
                     veclipos( j, k, nFix ) ] <- 1
               }
            }
         }
      } else {
         stop( "argument 'form' must be either 0 or 1" )
      }
   }
   return( MR )
}
