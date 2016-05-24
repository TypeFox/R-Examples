print.aidsMono <- function( x, header = TRUE, ... ) {
   if( header ) {
      cat( "\nChecking the monotonicity condition of an " )
      if( x$priceIndex == "TL" ) {
         cat( "Almost Ideal Demand System (AIDS):\n" )
      } else {
         cat( "Linear Approximate Almost Ideal Demand System (LA-AIDS):\n" )
      }
   }
   cat( "Monotonicity is fulfilled at " )
   cat( x$nMonoObs, "out of", x$nValidObs, "observations" )
   cat( " (", x$monoPercent, "%)\n", sep = "" )

   invisible( x )
}
