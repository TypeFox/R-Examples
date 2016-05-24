checkConsist.aidsEst <- function( object, observedShares = FALSE, ... ) {

   # price index for checking monotonicity
   if( object$method == "LA" ) {
      if( observedShares ) {
         priceIndex <- object$lnp
      } else {
         priceIndex <- object$priceIndex
      }
   } else if( object$method %in% c( "IL", "MK" ) ) {
      priceIndex <- "TL"
   } else {
      stop( "unknown element 'method' of argument 'object'" )
   }

   # shares for checking concavity
   if( observedShares ) {
      shareNames <- object$shareNames
   } else {
      shareNames <- NULL
   }


   aidsConsist( priceNames = object$priceNames,
      shareNames = shareNames,
      totExpName = object$totExpName,
      data = eval( object$call$data ),
      coef = object$coef,
      priceIndex = priceIndex,
      basePrices = object$basePrices,
      baseShares = object$baseShares )
}
