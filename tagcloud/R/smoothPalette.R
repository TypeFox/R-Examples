# turns a vector of numbers into a vector of colors


#' Replace a vector of numbers by a gradient of colors
#' 
#' Replace a vector of numbers by a vector of colors from a palette, such that
#' values correspond to the colors on a smooth gradient.
#' 
#' This function is used to map a continues numerical vector on an ordinal
#' character vector, in especially a vector of colors. Color palette can be
#' specified using an RColorBrewer palette name.
#' 
#' @param x A numeric vector
#' @param max Values of \code{x} larger than \code{max} will be replaced by
#' \code{max}
#' @param min Values of \code{x} smaller than \code{min} will be replaced by
#' \code{min}
#' @param pal Character vector containing the color gradient onto which the
#' numeric vector \code{x} will be mapped. By default, a gradient from white to
#' black is generated. If it is a single character value, it will be treated as
#' name of an RColorBrewer palette (see \code{\link{brewer.pal}}).
#' @param n Number of steps
#' @param palfunc Palette function returned by colorRampPalette
#' @param na.color NA values will be replaced by that color
#' @return A character vector of the same length as the numeric vector
#' \code{x}, containing the matching colors.
#' @author January Weiner <january.weiner@@gmail.com>
#' @seealso \code{\link{tagcloud}}
#' @keywords palette mapping
#' @examples
#' 
#' smoothPalette( 1:3 )
#' # will print:
#' # "#CCCCCC" "#666666" "#000000"
#' 
#' smoothPalette( 1:3, pal= "Blues" )
#' # will produce:
#' # "#F7FBFF" "#6BAED6" "#08306B"
#' 
#' x <- runif( 100 )
#' plot( 1:100, x, col= smoothPalette( x, pal= "BrBG" ), pch= 19 )
#' 
#' @export smoothPalette
smoothPalette <- function( x, pal= NULL, max= NULL, min= NULL, n= 9, palfunc= NULL,
                           na.color= "white" ) {

  if( ! missing( palfunc ) ) {
    pal <- palfunc( n )
  } else if( missing( pal ) ){
    pal <- colorRampPalette( c( "#cccccc", "black" ) )( n )
  } else {
    if( length( pal ) == 1 ) {
      pal <- try( brewer.pal( n, pal ), silent= TRUE )
      if( class( pal ) == "try-error" ) 
        stop( "palette is neither a vector nor a name of a RColorBrewer palette" )
    }
  }

  n <- length( pal )
  
  nas <- which( is.na( x ) )

  x[ nas ] <- mean( x, na.rm= T )

  if( missing( max ) ) max <- max( x )
  else x[ x >= max ] <- max

  if( missing( min ) ) min <- min( x )
  else x[ x <= min ] <- min

  ret <- findInterval( x, seq( min, max, length.out= n + 1 ), rightmost.closed= T )
  ret <- pal[ ret ]
  ret[ nas ] <- na.color
  return( ret )
}
