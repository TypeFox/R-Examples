## ===== calculation of elasticities from beta matrix ===
snqProfitEla <- function( beta, prices, quant, weights,
   scalingFactors = rep( 1, length( weights ) ),
   coefVcov = NULL, df = NULL ) {
   if( !is.matrix( beta ) ) {
      stop( "argument 'beta' must be a matrix" )
   }
   if( nrow( beta ) != ncol( beta ) ) {
      stop( "argument 'beta' must be a quadratic matrix" )
   }
   if( length( prices ) != length( quant ) ) {
      stop( "arguments 'prices' and 'quant' must have the same length" )
   }
   if( length( prices ) != length( weights ) ) {
      stop( "arguments 'prices' and 'weights' must have the same length" )
   }
   if( nrow( beta ) != length( prices ) ) {
      stop( "arguments 'prices' must have as many elements as",
         " argument 'beta' has rows" )
   }
   nNetput  <- ncol( beta )
   prices   <- unlist( prices ) * scalingFactors
   quant    <- unlist( quant ) / scalingFactors
   hessian  <- snqProfitHessian( beta, prices, weights )
   result   <- list()
   result$ela <- hessian * array( 1, c( nNetput ) ) %*% t( prices ) /
                  quant %*% t( array( 1, c( nNetput ) ) )

   quantNames <- .snqProfitQuantNames( quant, nNetput )
   priceNames <- .snqProfitPriceNames( prices, nNetput )

   dimnames( result$ela ) <- list( quantNames, priceNames )

   if( !is.null( coefVcov ) ) {
      jacobian <- snqProfitElaJacobian( beta, prices, quant, weights )
      betaIndex <- grep( "beta", rownames( coefVcov ) )
      betaVcov <- coefVcov[ betaIndex, betaIndex ]
      result$vcov <- jacobian %*% betaVcov %*% t( jacobian )
      result$stEr <- matrix( diag( result$vcov )^0.5, nrow = nNetput,
         byrow = TRUE )
      result$tval <- result$ela / result$stEr
      dimnames( result$stEr ) <- list( quantNames, priceNames )
      dimnames( result$tval ) <- list( quantNames, priceNames )
      if( !is.null( df ) ) {
         result$pval <- 2 * pt( abs( result$tval ), df,
            lower.tail = FALSE )
         dimnames( result$pval ) <- list( quantNames, priceNames )
      }
   }
   class( result ) <- "snqProfitEla"
   return( result )
}
