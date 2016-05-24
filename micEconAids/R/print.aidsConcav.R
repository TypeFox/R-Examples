print.aidsConcav <- function( x, header = TRUE, ... ) {
   if( header ) {
      cat( "\nChecking the concavity condition of an " )
      cat( "Almost Ideal Demand System (AIDS):\n" )
   }
   cat( "Concavity is fulfilled at " )
   cat( x$nConcavObs, "out of", x$nValidObs, "observations" )
   cat( " (", x$concavPercent, "%)\n", sep = "" )

   invisible( x )
}
