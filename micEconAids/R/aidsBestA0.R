aidsBestA0 <- function( priceNames, shareNames, totExpName,
      a0min = -50, a0max = 50, stoprange = 3, stopiter = 10,
      verbose = FALSE, ... ) {

   if( length( priceNames ) != length( shareNames ) ) {
      stop( "arguments 'priceNames' and 'shareNames' must have the same length" )
   }
   nGoods <- length( priceNames )

   if( a0min >= a0max ) stop( "a0min must be smaller than a0max" )

   deta0 <- function( a0, ... ) {
      estResult <- aidsEst( priceNames, shareNames, totExpName,
         method = "IL", alpha0 = a0, ... )
      det <- det( estResult$est$residCov )
      assign( "allValues", rbind( allValues, c( a0, det ) ),
         sys.frame( sys.parent( ) )  )
      if( verbose ) {
         cat( "a0:", a0, "-> det:", det, "(iterIL:", estResult$iterIL, ")\n" )
      }
      return( det )
   }

   a0 <- array( NA, c( 2, 4 ) )  # 1st row=alpha0; 2nd row=det(alpha0)
   a0[ 1, ] <- c( a0min, a0min + ( a0max - a0min ) / 3,
      a0max - ( a0max - a0min ) / 3, a0max )
   allValues <- NULL

   for( i in 1:4 ) {
      a0[2,i] <- deta0( a0[ 1, i ], ... )
   }
   iter <- 0
   while( ( which.min( a0[ 2, ] ) == 1 | which.min( a0[ 2, ]) == 4 ) &
      iter < stopiter ) {
      iter <- iter + 1
      if( which.min( a0[ 2, ] ) == 1 ) {
         a0[ , 2:4 ] <- a0[ , 1:3 ]
         a0[ 1, 1 ]  <- a0[ 1, 1 ] - ( a0max - a0min ) / 3
         a0[ 2, 1 ]  <- deta0( a0[ 1, 1 ], ... )
      } else {
         a0[ , 1:3 ] <- a0[ , 2:4 ]
         a0[ 1, 4 ]  <- a0[ 1, 4 ] + ( a0max - a0min ) / 3
         a0[ 2, 4 ]  <- deta0( a0[ 1, 4 ], ... )
      }
   }
   while( iter < stopiter & ( a0[ 1, 4 ] - a0[ 1, 1 ] ) > stoprange ) {
      iter <- iter + 1
      if( which.min( a0[ 2, ] ) == 2 ) {
         a0[ , 4 ] <- a0[ , 3 ]
         if( a0[ 1, 2 ] - a0[ 1, 1 ] >= a0[ 1, 3 ] - a0[ 1, 2 ]) {
            a0[ , 3 ]  <- a0[ , 2 ]
            a0[ 1, 2 ] <- ( a0[ 1, 1 ] + a0[ 1, 3 ] ) / 2
            a0[ 2, 2 ] <- deta0( a0[ 1, 2 ], ... )
         } else {
            a0[ 1, 3 ] <- ( a0[ 1, 2 ] + a0[ 1, 4 ] ) / 2
            a0[ 2, 3 ] <- deta0( a0[ 1, 3 ], ... )
         }
      } else if( which.min(a0[2,])==3) {
         a0[,1] <- a0[,2]
         if( a0[ 1, 4 ] - a0[ 1, 3 ] >= a0[ 1, 3 ] - a0[ 1, 2 ] ) {
            a0[ , 2 ]  <- a0[ , 3 ]
            a0[ 1, 3 ] <- ( a0[ 1, 2 ] + a0[ 1, 4 ] ) / 2
            a0[ 2, 3 ] <- deta0( a0[ 1, 3 ], ... )
         } else {
            a0[ 1, 2 ] <- ( a0[ 1, 1 ] + a0[ 1, 3 ] ) / 2
            a0[ 2, 2 ] <- deta0( a0[ 1, 2 ], ... )
         }
      } else {
         stop("minimum not between a0min and a0max")
      }
   }
   result <- list()
   result$alpha0 <- a0[ 1, which.min( a0[ 2, ] ) ]
   result$allValues <- allValues[ order( allValues[ , 1 ] ), ]
   colnames( result$allValues ) <- c( "a0", "det" )
   result$iter <- iter
   return( result )
}
