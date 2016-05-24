#' Draw a curved segment 
#'
#' Draws a curved segment from point \code{(x0,y0)} to \code{(x1,y1)}. The
#' segment is a framgent of a sinusoid, has a defined width and can either
#' have a single color or a color gradient.
#'
#'
#'
#'
#' @param x0 X coordinate of the starting point
#' @param y0 X coordinate of the starting point
#' @param x1 X coordinate of the end point
#' @param y1 X coordinate of the end point
#' @param width Width of the segment to plot
#' @param nsteps Number of polygons to use for the segments. The more, the
#' smoother the picture, but at the same time, the more time-consuming to
#' display.
#' @param col Color to use. Ignored if grad is not \code{NULL}.
#' @param grad Gradient to use. Can be anything that
#' \code{colorRampPalette} can understand.
#' @param lty Line type for drawing of the segment. Use \code{lty=0} for no line.
#' @param form "sin" for a sinusoidal segment. "line" for a straight segment.
#' @return no value is returned
#' @export
#' @examples 
#' # a DNA strand
#' plot.new()
#' par( usr= c( 0, 4, -2.5, 2.5 ) )
#'
#' w    <- 0.4
#' cols <- c( "blue", "green" )
#' init <- c( -0.8, -0.5 )
#' pos  <- c( 1, -1 )
#' step <- 0.5
#'
#' for( i in rep( rep( c( 1, 2 ), each= 2 ), 5 ) ) {
#'   curveseg( init[i], init[i] + step, pos[1], pos[2], width= w, col= cols[i] )
#'   init[i] <- init[i] + step
#'   pos <- pos * -1
#' }


curveseg <- function( x0, x1, y0, y1, width= 1, nsteps= 50, col= "#ffcc0066", grad= NULL, lty= 1, form= c( "sin", "line" ) ) {

  w <- width

  if( ! is.null( grad ) ) {
    grad <- colorRampPaletteAlpha( grad )( nsteps )
  } else {
    grad <- rep( col, nsteps )
  }

  form <- match.arg( form, c( "sin", "line" ) )

  if( form == "sin" ) {
    xx  <- seq( -pi/2, pi/2, length.out= nsteps )
    yy <- y0 + ( y1 - y0 ) * ( sin( xx ) + 1 ) / 2
    xx <- seq( x0, x1, length.out= nsteps )
  }

  if( form == "line" ) {
    xx <- seq( x0, x1, length.out= nsteps )
    yy <- seq( y0, y1, length.out= nsteps )
  }

  for( i in 1:(nsteps-1) ) {
    polygon( c( xx[i], xx[i+1], xx[i+1], xx[i] ),
             c( yy[i], yy[i+1], yy[i+1] + w, yy[i] + w ), col= grad[i], border= grad[i] )
             # c( yy[i], yy[i+1], yy[i+1] + w, yy[i] + w ), col= grad[i], lty= 0 )
    lines( c( xx[i], xx[i+1] ), c( yy[i], yy[i+1] ), lty= lty )
    lines( c( xx[i], xx[i+1] ), c( yy[i] + w, yy[i+1] + w ), lty= lty )
  }

}

