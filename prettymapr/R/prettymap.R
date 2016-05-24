#plot maps without space for axes


#' Plot A Pretty Map
#'
#' This function executes everything in \code{plotexpression}, then draws
#' north arrow and scale bar using \link{addnortharrow} and \link{addscalebar}.
#' Specify that plot is in a non lat/lon projection by passing \code{scale.plotepsg=...}
#' or \code{plotunit="m"}.
#'
#' @param plotexpression An expression to plot the map, can be in brackets. e.g.
#' \code{plot(stuff); text(places, "readme!")} or \code{{plot(stuff);
#' text(places, "readme!")}}
#' @param oma A vector of length 4 describing the outer margin area. See documentation
#' for \code{graphics::par}.
#' @param mai A vector of length 4 describing the margin area in inches. See documentation
#' for \code{graphics::par}.
#' @param drawbox \code{TRUE} if box should be drawn around map, \code{FALSE} otherwise.
#' @param box.lwd The line width of the box
#' @param drawscale \code{TRUE} if scalebar should be drawn, \code{FALSE} otherwise.
#' @param scale.plotunit The unit which the current plot is plotted in, one of \code{cm},
#' \code{m}, \code{km}, \code{in}, \code{ft}, \code{mi}. or \code{latlon}. This
#' parameter is optional if \code{plotepsg} is passed.
#' @param scale.plotepsg The projection of the current plot. If extents are valid lat/lons,
#' the projection is assumed to be lat/lon (EPSG:4326), or Spherical Mercator otherwise
#' (EPSG:3857). This is done to work seamlessly with OpenStreetMap packages.
#' @param scale.widthhint The fraction of the plottable width which the scale bar should
#' (mostly) occupy.
#' @param scale.unitcategory One of "metric" or "imperial"
#' @param scale.htin Height (in inches) of the desired scale bar
#' @param scale.padin A vector of length 2 determining the distance in inches between the scalebar
#' and the edge of the plottable area.
#' @param scale.style One of "bar" or "ticks".
#' @param scale.bar.cols If \code{style=="bar"}, the colors to be repeated to make the bar.
#' @param scale.lwd The line width to use when drawing the scalebar
#' @param scale.linecol The line color to use when drawing the scalebar
#' @param scale.tick.cex If \code{style=="ticks"}, the height of interior ticks.
#' @param scale.labelpadin The distance between the end of the scalebar and the label (inches)
#' @param scale.label.cex The font size of the label
#' @param scale.label.col The color of the label
#' @param scale.pos Where to align the scalebar. One of "bottomleft", "bottomright", "topleft",
#' or "topright".
#' @param drawarrow \code{TRUE} if north arrow should be drawn, \code{FALSE} otherwise
#' @param arrow.pos Where to align the north arrow. One of "bottomleft", "bottomright", "topleft",
#' or "topright".
#' @param arrow.padin A vector of length 2 determining the distance in inches between the scalebar
#' and the edge of the plottable area.
#' @param arrow.scale Scale the default north arrow to make it bigger or smaller
#' @param arrow.lwd The line width outlining the north arrow
#' @param arrow.border The line color outlining the north arrow
#' @param arrow.cols A vector of length 2 determining the two colors to be drawn for the north arrow
#' @param arrow.text.col Color of the "N"
#' @param title Plot title, or \code{NULL} if none is desired.
#' @param ... Further graphical parameters to set while executing plotting code
#'
#' @export
#'
#' @examples
#' prettymap(plot(1:5, 1:5, asp=1), scale.plotunit="cm", drawarrow=FALSE)
#' #add a title
#' prettymap(plot(1:5, 1:5, asp=1), title="My Plot")
#' \donttest{
#' library(maptools)
#' data(wrld_simpl)
#' prettymap({plot(wrld_simpl, xlim=c(-66.86, -59.75), ylim=c(43, 47.3))
#'            text(-62, 44, "Nova Scotia")
#'            text(-63, 47, "PEI")}, arrow.scale=1.1)
#'
#' #also works in non-lat/lon coordinate systems
#' prettymap(plot(1:1000, 1:1000, asp=1),
#'            scale.plotepsg=26920, drawarrow=FALSE) #specify plot is in UTM Zone 20N
#' }
#'
prettymap <- function(plotexpression, oma=c(0, 0, 0, 0),
                      mai=c(0, 0, 0, 0), drawbox=FALSE, box.lwd=1,
                      drawscale=TRUE, scale.pos="bottomleft", scale.htin=0.1,
                      scale.widthhint=0.25, scale.unitcategory="metric", scale.style="bar",
                      scale.bar.cols=c("black", "white"), scale.lwd=1, scale.linecol="black",
                      scale.padin=c(0.15, 0.15), scale.labelpadin=0.08, scale.label.cex=0.8,
                      scale.label.col="black", scale.plotunit=NULL, scale.plotepsg=NULL, scale.tick.cex=0.8,
                      drawarrow=FALSE, arrow.pos="topright", arrow.scale=1, arrow.padin=c(0.15, 0.15),
                      arrow.lwd=1, arrow.cols=c("white", "black"), arrow.border="black",
                      arrow.text.col="black", title=NULL, ...) {

  prevpars <- graphics::par(oma=oma, mai=mai, ...)
  tryCatch(expr={
    force(plotexpression)
    if(drawbox) graphics::box(lwd=box.lwd)
    if(drawscale) addscalebar(plotunit=scale.plotunit, pos=scale.pos, htin=scale.htin,
                           widthhint=scale.widthhint, unitcategory=scale.unitcategory, style=scale.style,
                           bar.cols=scale.bar.cols, lwd=scale.lwd, linecol=scale.linecol,
                           padin=scale.padin, labelpadin=scale.labelpadin, label.cex=scale.label.cex,
                           label.col=scale.label.col, plotepsg=scale.plotepsg, tick.cex=scale.tick.cex)
    if(drawarrow) addnortharrow(pos=arrow.pos, padin=arrow.padin, scale=arrow.scale, lwd=arrow.lwd,
                             cols=arrow.cols, border=arrow.border, text.col=arrow.text.col)
    }, error=function(e) {
    message(paste("Error occured while plotting: ", e))
  }, finally={graphics::par(prevpars)})

  if(!is.null(title)) {
    prevpars <- graphics::par(mar=c(0,0,1,0))
    graphics::title(title)
    graphics::par(prevpars)
  }
}
