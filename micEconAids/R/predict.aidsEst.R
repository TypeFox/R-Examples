predict.aidsEst <- function( object, newdata = NULL,
      observedShares = FALSE, ... ) {

   if( is.null( newdata ) ) {
      newdata <- eval( object$call$data )
   }

   # price index for checking monotonicity
   if( object$method == "LA" ) {
      if( observedShares ) {
         priceIndex <- aidsPx( priceIndex = object$priceIndex,
            priceNames = object$priceNames,
            data = newdata,
            shareNames = object$shareNames,
            base = list( prices = object$basePrices, shares = object$baseShares ),
            coef = coef( object ), shifterNames = object$shifterNames )
      } else {
         priceIndex <- object$priceIndex
      }
   } else if( object$method %in% c( "IL", "MK" ) ) {
      priceIndex <- "TL"
   } else {
      stop( "unknown element 'method' of argument 'object'" )
   }

   result <- aidsCalc( priceNames = object$priceNames,
      totExpName = object$totExpName,
      coef = coef( object ),
      data = newdata,
      priceIndex = priceIndex,
      basePrices = object$basePrices,
      baseShares = object$baseShares )

   return( result )
}