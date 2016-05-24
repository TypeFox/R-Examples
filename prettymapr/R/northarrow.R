
#' Plot North Arrow
#'
#' Plot a north arrow (pointing directly "up") positioned based on current
#' plot extents.
#'
#' @param pos Where to align the north arrow. One of "bottomleft", "bottomright", "topleft",
#' or "topright".
#' @param padin A vector of length 2 determining the distance in inches between the scalebar
#' and the edge of the plottable area.
#' @param scale Scale the default north arrow to make it bigger or smaller
#' @param lwd The line width outlining the north arrow
#' @param border The line color outlining the north arrow
#' @param cols A vector of length 2 determining the two colors to be drawn for the north arrow
#' @param text.col Color of the "N"
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(maptools)
#' data(wrld_simpl)
#' plot(wrld_simpl)
#' addnortharrow()
#' }
#' plot(1:5, 1:5, asp=1)
#' addnortharrow()
#'
addnortharrow <- function(pos="topright", padin=c(0.15, 0.15), scale=1, lwd=1, border="black",
                       cols=c("white", "black"), text.col="black") {
  extents <- graphics::par("usr")

  arrow1.x <- c(-0.5, 0, 0)
  arrow1.y <- c(-0.5, 0.5, 0)

  arrow2.x <- -1 * arrow1.x
  arrow2.y <- arrow1.y

  htin <- 0.6 * scale
  wdin <- 0.6 * scale
  text.cex <- 1.25 * scale

  bottomin <- graphics::grconvertY(extents[3], from="user", to="inches")
  leftin <- graphics::grconvertX(extents[1], from="user", to="inches")
  topin <- graphics::grconvertY(extents[4], from="user", to="inches")
  rightin <- graphics::grconvertX(extents[2], from="user", to="inches")
  textheight <- graphics::strheight("N", units="user", cex=text.cex)

  if(pos=="bottomleft") {
    x <- graphics::grconvertX(leftin+padin[1], from="inches", to="user")
    y <- graphics::grconvertY(bottomin+padin[2], from="inches", to="user") + textheight/2
    adj <- c(0,0)
  } else if(pos=="topleft") {
    x <- graphics::grconvertX(leftin+padin[1], from="inches", to="user")
    y <- graphics::grconvertY(topin-padin[2], from="inches", to="user")
    adj <- c(0,1)
  } else if(pos=="topright") {
    x <- graphics::grconvertX(rightin-padin[1], from="inches", to="user")
    y <- graphics::grconvertY(topin-padin[2], from="inches", to="user")
    adj <- c(1,1)
  } else if(pos=="bottomright") {
    x <- graphics::grconvertX(rightin-padin[1], from="inches", to="user")
    y <- graphics::grconvertY(bottomin+padin[2], from="inches", to="user") + textheight/2
    adj <- c(1,0)
  }

  yin <- graphics::grconvertY(y, from="user", to="inches")
  xin <- graphics::grconvertX(x, from="user", to="inches")
  ht <- graphics::grconvertY(yin+htin, from="inches", to="user") - y
  wd <- graphics::grconvertX(xin+wdin, from="inches", to="user") - x

  graphics::polygon(x-adj[1]*wd+arrow1.x*wd+wd/2, y-adj[2]*ht+arrow1.y*ht+ht/2,
                   border=border, lwd=lwd, col=cols[1])
  graphics::polygon(x-adj[1]*wd+arrow2.x*wd+wd/2, y-adj[2]*ht+arrow2.y*ht+ht/2,
                   border=border, lwd=lwd, col=cols[2])
  graphics::text(x-adj[1]*wd+wd/2, y-adj[2]*ht+ht*0, "N",
                 cex=text.cex, adj=c(0.5, 0.5), col=text.col)

}
