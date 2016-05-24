elas.aidsEst <- function( object, method = NULL, observedShares = FALSE, ... ) {

   # specify default value for argument method
   if( is.null ( method ) ) {
      if( substr( object$method, 1, 2 ) == "LA" ) {
         method <- "Ch"
      } else if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) ) {
         method <- "AIDS"
      }
   }

   # test reasonability of argument method
   if( substr( object$method, 1, 2 ) %in% c( "MK", "IL" ) &&
         method != "AIDS" ) {
      warning( paste( "It does not make sense to calculate the elasticities",
         " of a (non-linear) AIDS model with method '", method, "'",
         sep = "" ) )
   }

   if( object$method == "LA" ) {
      priceIndex <- object$priceIndex
   } else if( object$method %in% c( "IL", "MK" ) ) {
      priceIndex <- "TL"
   } else {
      stop( "unknown element 'method' of argument 'object'" )
   }

   if( observedShares ) {
      shares <- object$wMeans
      totExp <- NULL
   } else {
      shares <- NULL
      totExp <- object$xtMean
   }

   # calculate demand elasticities
   result  <- aidsElas( coef = coef( object ),
      shares = shares, prices = object$pMeans, totExp = totExp,
      method = method,
      priceIndex = priceIndex,
      basePrices = object$basePrices,
      baseShares = object$baseShares,
      priceNames = object$priceNames,
      coefCov = vcov( object ), df = df.residual( object ), ... )

   return( result )
}

