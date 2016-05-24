createSystemfitModel <- function( nEq, nRegEq, nObs, coef = NULL, sigma = NULL ){

   result <- list()
   nCoef <- nEq * ( nRegEq + 1 )

   if( is.null( coef ) ){
      coef <- round( rnorm( nCoef ), 1 )
   } else {
      if( length( coef ) != nCoef ){
         stop( "argument 'coef' must have length ", nCoef )
      }
   }

   if( is.null( sigma ) ){
      sigma <- matrix( rnorm( nEq^2 ), nEq, nEq )
      sigma <- crossprod( sigma )
   } else {
      if( all.equal( dim( sigma ), c( nEq, nEq ) ) != TRUE ){
         stop( "argument 'sigma' must be a ", nEq, "x", nEq, " matrix" )
      }
   }


   disturbances <- mvrnorm( nObs, rep( 0, nEq ), sigma )
   result$data <- data.frame( obsNo = c( 1:nObs ) )
   result$formula <- list()

   for( eqNo in 1:nEq ){
      result$formula[[ eqNo ]] <- as.formula( paste( "y.", eqNo, " ~ ",
         paste( paste( "x.", eqNo, ".", sep = "" ),
         c( 1:nRegEq ), sep = "", collapse = " + " ), sep = "" ) )
      for( regNo in 1:nRegEq ){
         result$data[[ paste( "x", eqNo, regNo, sep = "." ) ]] <-
            rnorm( nObs )
      }
      result$data[[ paste( "y", eqNo, sep = "." ) ]] <- 0
      xMatEq <- model.matrix( result$formula[[ eqNo ]], result$data )
      result$data[[ paste( "y", eqNo, sep = "." ) ]] <-
         drop( xMatEq %*% coef[ ( ( eqNo - 1 ) * ( nRegEq + 1 ) + 1 ):
            ( eqNo * ( nRegEq + 1 ) ) ] ) + disturbances[ , eqNo ]
   }

   result$coef <- coef
   result$sigma <- sigma

   return( result )
}