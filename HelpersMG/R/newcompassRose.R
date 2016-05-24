#' newcompassRose Display a compass rose
#' @title Display a compass rose
#' @author modified from Jim Lemon; See compassRose {sp}
#' @return none
#' @param x The position of the center of the compass rose in user units.
#' @param y The position of the center of the compass rose in user units.
#' @param rot Rotation for the compass rose in degrees. See Details.
#' @param cex The character expansion to use in the display.
#' @param col The color of text
#' @param col.arrows.light The color of lighter lines
#' @param col.arrows.dark The color of darker lines
#' @description Displays a basic compass rose, usually to orient a map.\cr
#' newcompassRose displays a conventional compass rose at the position requested.\cr 
#' The size of the compass rose is determined by the character expansion, 
#' as the central "rose" is calculated relative to the character size.\cr
#' Rotation is in degrees counterclockwise.
#' @examples
#' \dontrun{
#' library(HelpersMG)
#' require("maps")
#' map("world", "China")
#' newcompassRose(x=110, y=35, col.arrows.light="grey")
#' }
#' @export


newcompassRose <- function (x, y, rot = 0, cex = 1, col="black", 
	col.arrows.light="white", col.arrows.dark="black") 
{
  oldcex <- par(cex = cex)
  mheight <- strheight("M")
  xylim <- par("usr")
  plotdim <- par("pin")
  xmult <- (xylim[2] - xylim[1])/(xylim[4] - xylim[3]) * plotdim[2]/plotdim[1]
  point.angles <- seq(0, 2 * pi, by = pi/4) + pi * rot/180
  crspans <- rep(c(mheight * 3, mheight/2), length.out = 9)
  xpoints <- cos(point.angles) * crspans * xmult + x
  ypoints <- sin(point.angles) * crspans + y
  for (point in 1:8) {
    pcol <- ifelse(point%%2, col.arrows.dark, col.arrows.light)
    polygon(c(xpoints[c(point, point + 1)], x), c(ypoints[c(point, 
                                                            point + 1)], y), col = pcol)
  }
  txtxpoints <- cos(point.angles[c(1, 3, 5, 7)]) * 1.2 * crspans[1] * 
    xmult + x
  txtypoints <- sin(point.angles[c(1, 3, 5, 7)]) * 1.2 * crspans[1] + 
    y
  text(txtxpoints, txtypoints, c("E", "N", "W", "S"), col=col)
  par(oldcex)
}
