print.summary.front41Output <- function( x, efficiencies = FALSE, ... ) {
   cat( "\nStochastic Frontier Analysis with FRONTIER 4.1\n" )
   cat( "Model type:", x$modelTypeName, "\n" )
   cat( "Function type:", x$functionTypeName, "\n" )
   cat( "\nML Estimates:\n" )
   printCoefmat( coef( x ) )
   if( efficiencies ) {
      cat( "\nEfficiency Estimates:\n" )
      print( x$efficiency )
   }
   cat( "\nMean Efficiency:", x$meanEfficiency, "\n" )
   invisible( x )
}
