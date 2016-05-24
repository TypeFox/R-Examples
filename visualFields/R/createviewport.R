createviewport <- function( name, left, top, width, height, pwidth = 8.27,
                            pheight = 11.69 ) {

  left   <- left / pwidth
  width  <- width / pwidth
  top    <- ( pheight - top ) / pheight
  height <- height / pheight

  return( viewport( x = left, y = top, width = width, height = height, just = c( "left", "top" ), name = name, clip = FALSE ) )

}