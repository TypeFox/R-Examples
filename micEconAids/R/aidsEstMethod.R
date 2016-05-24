.aidsEstMethod <- function( method, priceIndex ) {

   if( priceIndex == "S" ) {
      result <- "Stone Index"
   } else if( priceIndex == "SL" ) {
      result <- "lagged Stone Index"
   } else if( priceIndex == "P" ) {
      result <- "Paasche Index"
   } else if( priceIndex == "L" ) {
      result <- "Laspeyres Index"
   } else if( priceIndex == "Ls" ) {
      result <- "simplified Laspeyres Index"
   } else if( priceIndex == "T" ) {
      result <- "Tornqvist Index"
   } else {
      result <- "unknown price index"
   }
   if( substr( method, 1, 2 ) == "LA" ) {
      result <- paste( "Linear Approximation (LA) with ",
         result, " (", priceIndex, ")\n", sep = "" )
   } else if( substr( method, 1, 2 ) %in% c( "MK", "IL" ) ) {
      result <- paste( "'Iterated Linear Least Squares Estimator' (IL)\n",
         "(starting with ", result, ", ", priceIndex, ")\n", sep = "" )
   } else {
      result <- paste( "unknown method with ", result, " (",
         priceIndex, ")\n", sep = "" )
   }

   return( result )
}
