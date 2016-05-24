#' Levelplot of Raw Data
#' 
#' Plots raw data values without topography (height adjustment).
#' The raw data values have not been inverted yet.
#' @param Profile profile.
#' @param xlab label for x-axes.
#' @param ylab label for y-axes.
#' @param main title to be plotted.
#' @param trafo function to transform raw data values (default: log).
#' @param col vector of colors.
#' @param ... lattice levelplot arguments.
#' @export
#' @examples 
#' data(sinkhole)
#' 
#' levelplotRaw(sinkhole@profiles[[1]])
#' levelplotLegendLabel()
#' 
#' levelplotRaw(sinkhole@profiles[[2]])
#' levelplotLegendLabel()
#' 
#' levelplotRaw(sinkhole@profiles[[3]])
#' levelplotLegendLabel()
levelplotRaw <- function(Profile, xlab="Length [m]", ylab="Depth [m]",
                         main=paste(Profile@title, "without topography (raw data)"), 
                         col=colors, trafo=log, ...) {
  levelplot(trafo(Profile@rawData@seaLevel$val) ~ Profile@rawData@seaLevel$dist * (-1*Profile@rawData@seaLevel$depth), 
            col.regions = colorRampPalette(col), 
            xlab=xlab, ylab=ylab, main=main, ...)
}