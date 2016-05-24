compPlot <- function( x, y, lim = NULL, ... ) {

   xyRange <- range( x, y, na.rm = TRUE, finite = TRUE )

   if( is.null( lim ) ) {
      lim <- xyRange
   } else {
      if( length( lim ) != 2 ) {
         stop( "argument 'lim' must be a vector of two elements" )
      }
      if( is.na( lim[1] ) ) {
         lim[1] <- xyRange[1]
      }
      if( is.na( lim[2] ) ) {
         lim[2] <- xyRange[2]
      }
      if( lim[1] >= lim[2] ) {
         stop( "the first element of argument 'lim' must be smaller",
            " than the second element" )
      }
      if( lim[1] > xyRange[1] |  lim[2] < xyRange[2] ) {
         warning( "some data points are outside the print area" )
      }
   }
   
   # code taken from plot.default()
   xlabel <- deparse(substitute(x))
   ylabel <- deparse(substitute(y))

   argList <- list( ... )
   
   argList$x <- x
   argList$y <- y
   argList$xlim <- lim
   argList$ylim <- lim
   if( ! "xlab" %in% names (argList) ) {
      argList$xlab <- xlabel 
   }
   if( ! "ylab" %in% names (argList) ) {
      argList$ylab <- ylabel 
   }
   
   do.call( plot.default, argList )

   abline( 0, 1 )

   invisible( xyRange )
}
