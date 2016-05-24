aidsMono <- function( priceNames, totExpName, coef, data,
      priceIndex = "TL", basePrices = NULL, baseShares = NULL ) {

   if( is.character( priceIndex ) ) {
      if( is.null( coef$alpha0 ) && priceIndex == "TL" ) {
         stop( "argument 'coef' must have element 'alpha0'" )
      }
   }

   result <- list()
   nGoods <- length( priceNames )
   nObs <- nrow( data )

   # fitted shares
   fitted <- aidsCalc( priceNames = priceNames, totExpName = totExpName,
      coef = coef, data = data, priceIndex = priceIndex,
      basePrices = basePrices, baseShares = baseShares )


   # testing for monotonicity
   result$monotony <- rep( TRUE, nObs )
   if( !is.null( row.names( data ) ) ) {
      names( result$monotony ) <- row.names( data )
   }
   for( t in 1:nObs ) {
      result$monotony[ t ] <- ( min( fitted$shares[ t, ] ) >= 0 )
   }

   result$nValidObs <- sum( !is.na( result$monotony ) )
   result$nMonoObs <- sum( result$monotony, na.rm = TRUE )
   result$monoPercent <- 100 * result$nMonoObs / result$nValidObs
   if( is.character( priceIndex ) ) {
      result$priceIndex <- priceIndex
   } else {
      result$priceIndex <- "numeric"
   }

   class( result ) <- "aidsMono"
   return( result )
}
