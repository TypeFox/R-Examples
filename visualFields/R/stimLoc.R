stimLoc <- function( perimetry, pattern, eye, txtfont = "mono", pointsize = 7,
                     xminmax = 29, yminmax = 29 ) {
  evaltxt <- paste( perimetry, "locmap$", pattern, sep = "" )
  locmap  <- eval( parse( text = evaltxt ) )
# left/right eye
  if( eye == "OS" ) {
    locmap$xod <- -locmap$xod
  }
  if( perimetry == "fdp" ) {
    stiSymbol <- "square"
    outerDimensions <- locmap$size
  } else {
    stiSymbol <- "circle"
    outerDimensions <- locmap$size / 2
  }
  plot( locmap$xod, locmap$yod, type = "n", axes = FALSE, xlab = "", ylab = "", asp = 1, xlim = c( -xminmax, xminmax ), ylim = c( -yminmax, yminmax ) )
  axis( 1, pos = 0, labels = FALSE, lwd = 0.5, lwd.ticks = 0 )
  axis( 2, pos = 0, las = 1, lwd = 0.5, lwd.ticks = 0, at = c( -yminmax, yminmax ), labels = c( yminmax, yminmax ), hadj = -0.5, padj = c( 2, -1 ) )

  evaltxt <- paste( "symbols( locmap$xod, locmap$yod, " , stiSymbol, " = outerDimensions, add = TRUE, lwd = 0.5, inches = FALSE )", sep = "" )
  eval( parse( text = evaltxt ) )
}
