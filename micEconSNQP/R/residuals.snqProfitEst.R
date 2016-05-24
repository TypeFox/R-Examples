residuals.snqProfitEst <- function( object, scaled = TRUE, ... ) {

   nNetput <- length( object$pMeans )
   nFixed  <- length( object$fMeans )
   nObs    <- nrow( object$data )

   result <- data.frame( profit0 = rep( 0, nObs ) )
   for( i in 1:nNetput ) {
      result[[ object$quantNames[ i ] ]] <-
         object$data[[ object$quantNames[ i ] ]] /
            object$scalingFactors[ i ]^( scaled ) -
         object$fitted[[ object$quantNames[ i ] ]] *
            object$scalingFactors[ i ]^( !scaled )
      result$profit0 <- result$profit0 +
         object$data[[ object$quantNames[ i ] ]] *
         object$data[[ object$priceNames[ i ] ]]
   }
   result$profit <- result$profit0 - object$fitted$profit
   result$profit0 <- NULL

   return( result )
}

## the same for snqProfitImposeConvexity
residuals.snqProfitImposeConvexity <- function( object, scaled = TRUE, ... ) {

   result <- residuals.snqProfitEst( object, scaled = scaled, ... )

   return( result )
}
