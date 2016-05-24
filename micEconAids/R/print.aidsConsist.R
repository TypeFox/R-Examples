print.aidsConsist <- function( x, ... ) {
   cat( "\nChecking theoretical consistency of an " )
   if( x$mono$priceIndex == "TL" ) {
      cat( "Almost Ideal Demand System (AIDS):\n" )
   } else {
      cat( "Linear Approximate Almost Ideal Demand System (LA-AIDS):\n" )
   }

   # Addinp-up
   cat( "The adding-up condition is" )
   if( !x$addingUp ) {
      cat( " NOT" )
   }
   cat( " fulfilled\n" )

   # homogeneity
   cat( "The homogeneity condition is" )
   if( !x$homogeneity ) {
      cat( " NOT" )
   }
   cat( " fulfilled\n" )

   # symmetry
   cat( "The symmetry condition is" )
   if( !x$symmetry ) {
      cat( " NOT" )
   }
   cat( " fulfilled\n" )

   # monotonicity
   print( x$mono, header = FALSE )

   # concavity
   if( !is.null( x$concav ) ) {
      print( x$concav, header = FALSE )
   }
   invisible( x )
}
