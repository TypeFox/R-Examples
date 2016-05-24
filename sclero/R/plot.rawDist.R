#' @title Plot rawDist object
#' @description Plots a map of \code{\link[=convert.ijdata]{rawDist}} object.
#' 
#' @param x \code{\link[=convert.ijdata]{rawDist}} object
#' @param ... Arguments to be passed to other methods, such as \link[=par]{graphical parameters}.
#' @param sample.name A character argument specifying the sample name to be plotted as an overall title for the plot (\code{\link[=plot.default]{main}}). Defaults to \code{"keep"} meaning that the sample name will be extracted from the \code{\link[=convert.ijdata]{rawDist}} object. The plot title can be omitted by specifying \code{sample.name = NULL}.
#' @param spot.type A character argument with three possible levels (\code{"id"}, \code{"value"}, and \code{"idvalue"}) indicating how sample spots should be plotted. Defaults to \code{"id"}, which plots sample spot numbers within open circles. The size of the circles can be controlled using the \code{spot.size} argument. The option \code{"value"} results to a sample map where the color of circles is related to a value through \code{\link{assign.value}} function. The color scale can be set using the \code{color.palette} argument, and size of the symbols (pch = 21) and through the \code{spot.size} argument. The option \code{"idvalue"} combines \code{"id"} and \code{"value"} leading to a sample map with sample spot numbers.
#' @param spot.size An integer or a character argument with value \code{"actual"} indicating the size (\code{\link[=par]{cex}}) of points. If \code{"actual"}, the actual size and shape of sample spots will be plotted. In this case, \link[=assign.size]{sample spot size information} is required. Defaults to 2 meaning that sample spots are plotted as points with pch = 21 and cex = 2.
#' @param spot.color A vector with equal length to number of spot sequences defining the color of sample spot labels. If \code{NULL} (default) a preset set of colors will be used.
#' @param main.type A character argument with four possible levels (\code{"all"}, \code{"axis"}, \code{"ends"}, and \code{"none"}) indicating how the distance / main axis should be plotted. Defaults to \code{"all"} indicating that both the main axis and end points should be plotted. If \code{"axis"} only the main axis will be plotted. If \code{"ends"} only the end points will be plotted, and if \code{"none"} the main axis intormation will not be plotted.
#' @param color.palette color palette used for "value" and "idvalue" \code{spot.type} options. Passed to \code{\link[grDevices]{colorRampPalette}}. 
#' @param highlight.gbs A character vector specifying the names of growth bands to be highlighted (i.e. colored with a different color than "darkgrey"). If \code{NULL} (default) all growth bands will be drawn using the standard color.
#' @param highlight.col A character argument specifying the color to be used in growth band highlighting (\code{highlight.gbs}).
#' @details The \pkg{sclero} package currently uses the \pkg{graphics} package distributed with R for plotting (see \code{\link[graphics]{plot}}). Plotting sample maps is carried out by the \code{sclero:::samplemap} function, which works as an internal function and therefore has not been exported. Users willing to modify \pkg{sclero} plots beyond the flexibility allowed by the \code{\link{plot.rawDist}} function are instructed to modify the \code{\link{samplemap}} function, which consists of standard R graphics syntax. It should be noted that the \code{\link{samplemap}} (and therefore also the \code{\link{plot.rawDist}}) function calls for the {\code{\link[graphics]{layout}}} function every time the arguments \code{spot.type = "value"} or \code{spot.type = "idvalue"} are used. Consequently, the graphics window is divided into two regions that might cause issues when combining \pkg{sclero} plots with other graphics. The users are adviced to consider the graphics window resetting procedure specified in {\code{\link[graphics]{layout}}} examples.
#' 
#' Because the function plots a sample map, the \strong{aspect ratio} is forced to 1 and cannot be changed. If this causes troubles when trying to set the axis limits (\code{ylim} and \code{xlim}), try resizing the graphics window.
#' @author Mikko Vihtakari
#' @seealso \code{\link{convert.ijdata}} for converting the coordinate information to \link[spatstat]{spatstat} point patterns. 
#' 
#' \code{\link{read.ijdata}} for reading zip files containing ImageJ ROIs.
#' 
#' \code{\link{plot.spotDist}} for plotting \code{\link[=spot.dist]{spotDist}} objects.
#' 
#' \code{\link[=plot]{plot.default}} and other methods; \code{\link[=points]{points}}, \code{\link[=lines]{lines}}, \code{\link[=par]{par}}.
#' 
#' @examples data(shellspots)
#' shell_map <- convert.ijdata(shellspots)
#' plot(shell_map)
#' 
#' @import spatstat 
#' @importFrom grDevices colorRampPalette
#' @importFrom stats na.omit
#' @importFrom graphics axis layout par plot points text
#' @export

plot.rawDist <- function(x, ..., sample.name = "keep", spot.type = "id", spot.size = 2, spot.color = NULL, main.type = "all", color.palette = colorRampPalette(c("blue", "cyan", "yellow", "red"), bias=1)(100), highlight.gbs = NULL, highlight.col = "red"){

## Plot the sample map
output <- samplemap(x = x, ..., sname = sample.name, sptype = spot.type, size = spot.size, scol = spot.color, mtype = main.type, colpalette = color.palette, hlight = highlight.gbs, hlcol = highlight.col)

## Add legend
  
if(spot.type == "value" | spot.type == "idvalue") {
par(mar=c(1,0.5,3,2))
   plot(output$colmap, vertical = TRUE, las = 2, main = x$value.name, cex = 0.1)}
}