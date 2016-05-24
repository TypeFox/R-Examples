checkNames <- function( testNames, allNames ) {
   inAllNames   <- testNames %in% allNames
   if( !all( inAllNames ) ) {
      stop( "object(s) '",
         paste( testNames[ !inAllNames ], collapse = "', '" ),
         "' not found." )
   }
}

