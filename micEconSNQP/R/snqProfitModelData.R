.snqProfitModelData <- function( data, weights, priceNames, quantNames, fixNames,
   instNames, form, netputScale = rep( 1, length( weights ) ),
   fixedScale = rep( 1, length( fixNames ) ) ){

   nNetput <- length( quantNames )  # number of netputs
   nFix    <- length( fixNames )  # number of fixed inputs
   nIV     <- length( instNames )  # number of fixed inputs
   nObs    <- nrow( data )      # number of observations

   ## price index for normalization
   result <- data.frame( nr = 1:nObs, normPrice = 0 )
   for( i in 1:nNetput ) {
      result$normPrice <- result$normPrice +
         data[[ priceNames[ i ] ]] * netputScale[ i ] * weights[ i ]
   }

   ## real/normalized netput prices and netput quantities
   for( i in 1:nNetput ) {
      result[[ paste( "pr", as.character( i ), sep = "" ) ]] <-
         data[[ priceNames[ i ] ]] * netputScale[ i ] / result$normPrice
      result[[ paste( "q", as.character( i ), sep = "" ) ]] <-
         data[[ quantNames[ i ] ]] / netputScale[ i ]
   }

   ## quadratic netput prices
   for( i in 1:nNetput ) {
      for( j in 1:nNetput ) {
         for( k in 1:nNetput ) {
            result[[ paste( "pq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
               -0.5 * weights[ i ] * data[[ priceNames[ j ] ]] * netputScale[ j ] *
               data[[ priceNames[ k ] ]] * netputScale[ k ] / result$normPrice^2
         }
      }
   }

   ## quasi-fix inputs
   if( nFix > 0 ) {
      for( i in 1:nFix ) {
         result[[ paste( "f", as.character( i ), sep = "" ) ]] <-
            data[[ fixNames[ i ] ]] / fixedScale[ i ]
      }
      ## quadratic quasi-fix inputs
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            for( k in 1:nFix ) {
               result[[ paste( "fq", as.character( i ), ".", as.character( j ), ".",
               as.character( k ), sep = "" ) ]] <-
                  0.5 * ifelse( form == 0, weights[ i ], 1 ) *
                  ( data[[ fixNames[ j ] ]] / fixedScale[ j ] ) *
                  ( data[[ fixNames[ k ] ]] / fixedScale[ k ] )
            }
         }
      }
   }

   ## instrumental variables
   if( nIV > 0 ) {
      for( i in 1:nIV ) {
         result[[ paste( "iv", as.character( i ), sep = "" ) ]] <-
            data[[ instNames[ i ] ]] / mean( data[[ instNames[ i ] ]] )
      }
   }

   ## model data for the profit function (not for the netput equations)
   ## current netput prices
   for( i in 1:nNetput ) {
      result[[ paste( "tp", as.character( i ), sep = "" ) ]] <-
         data[[ priceNames[ i ] ]] * netputScale[ i ]
   }
   ## quadratic netput prices
   for( i in 1:nNetput ) {
      for( j in 1:nNetput ) {
         result[[ paste( "tpq", as.character( i ), ".", as.character( j ),
            sep = "" ) ]] <- 0.5 * data[[ priceNames[ i ] ]] * netputScale[ i ] *
            data[[ priceNames[ j ] ]] * netputScale[ j ] / result$normPrice
      }
   }
   if( nFix > 0 ) {
      ## netput prices x quasi-fix inputs
      for( i in 1:nNetput ) {
         for( j in 1:nFix ) {
            result[[ paste( "tpf", as.character( i ), ".", as.character( j ),
               sep = "" ) ]] <- data[[ priceNames[ i ] ]] * netputScale[ i ] *
               data[[ fixNames[ j ] ]] / fixedScale[ j ]
         }
      }
      ## quadratic quasi-fix inputs
      if( form == 0 ) {
         for( i in 1:nFix ) {
            for( j in 1:nFix ) {
               result[[ paste( "tfq", as.character( i ), ".", as.character( j ),
                  sep = "" ) ]] <- 0.5 * result$normPrice *
                  ( data[[ fixNames[ i ] ]] / fixedScale[ i ] ) *
                  ( data[[ fixNames[ j ] ]] / fixedScale[ j ] )
            }
         }
      } else {
         for( i in 1:nNetput ) {
            for( j in 1:nFix ) {
               for( k in 1:nFix ) {
                  result[[ paste( "tfq", as.character( i ), ".",
                     as.character( j ), ".", as.character( k ), sep = "" ) ]] <-
                     0.5 * data[[ priceNames[ i ] ]] * netputScale[ i ] *
                     ( data[[ fixNames[ j ] ]] / fixedScale[ j ] ) *
                     ( data[[ fixNames[ k ] ]] / fixedScale[ k ] )
               }
            }
         }
      }
   }
   result$profit <- 0
   for( i in 1:nNetput ) {
      result$profit <- result$profit +
         data[[ priceNames[ i ] ]] * data[[ quantNames[ i ] ]]
   }

   return( result )
}
