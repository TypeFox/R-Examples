library(maptools)
library(ggplot2)

my.shp <- readShapePoly("/media/My Passport/data/GIS/USA/USA_State_2000/fe_2007_us_state00")

str(state_map)


#' Turns a SPDF into a data.frame that ggplot2 handles nicely
#' @param model Ignored
#' @param dataset The SpatialPolygonsDataFrame to parse
#' @param include.data.frame Whether to include the data component of the SPDF (defaults to yes, but for large polygons this may take up too much memory)
#' @return A data.frame suitable for use with ggplot2
fortify.SpatialPolygonsDataFrame <- function(model, dataset, include.data.frame=TRUE, ...) {
  
}