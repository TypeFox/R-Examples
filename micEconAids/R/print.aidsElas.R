print.aidsElas <- function( x, ... ) {

   cat( "\nDemand Elasticities " )

   if( x$priceIndex == "TL" ) {
      pxName <- "translog"
   } else if( x$priceIndex %in% c( "S", "SL" ) ) {
      pxName <- "Stone"
   } else if( x$priceIndex == "P" ) {
      pxName <- "Paasche"
   } else if( x$priceIndex %in% c( "L", "Ls" ) ) {
      pxName <- "Laspeyres"
   } else if( x$priceIndex == "T" ) {
      pxName <- "Tornqvist"
   } else {
      pxName <- "unknown"
   }

   if( x$method %in% c( "GA", "B1" ) ) {
      methodName <- "formulas of Green and Alston / Buse (first set)"
   } else if( x$method %in% c( "B2" ) ) {
      methodName <- "formulas of Buse (second set)"
   } else if( x$method %in% c( "Go", "Ch" ) ) {
      methodName <- "formulas of Goddard / Chalfant"
   } else if( x$method == "EU" ) {
      methodName <- "formulas of Eales and Unnevehr"
   } else if( x$method == "AIDS" ) {
      methodName <- "original AIDS formulas"
   } else {
      methodName <- paste( "unknown formula '", x$method, "'", sep = "" )
   }

   if( x$method %in% c( "AIDS", "EU" ) ) {
      cat( "(", methodName, ")\n", sep = "" )
   } else {
      cat( "(", methodName, " for ", pxName, " price index)\n", sep = "" )
   }

   cat( "Expenditure Elasticities\n" )
   print( x$exp )
   cat( "\nMarshallian (uncompensated) Price Elasticities\n" )
   print( x$marshall )
   cat( "\nHicksian (compensated) Price Elasticities\n" )
   print( x$hicks )
   invisible( x )
}
