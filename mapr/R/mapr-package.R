#' Vizualize species occurrence data
#'
#' @section Many inputs:
#' All functions take the following kinds of inputs:
#' \itemize{
#'  \item An object of class \code{occdat}, from the package \pkg{spocc}. An object of
#'  this class is composed of many objects of class \code{occdatind}
#'  \item An object of class \code{occdatind}, from the package \pkg{spocc}
#'  \item An object of class \code{gbif}, from the package \pkg{rgbif}
#'  \item An object of class \code{data.frame}. This data.frame can have any columns,
#'  but must include a column for taxonomic names (e.g., \code{name}), and for latitude
#'  and longitude (we guess your lat/long columns, starting with the default
#'  \code{latitude} and \code{longitude})
#'  \item An object of class \code{SpatialPoints}
#'  \item An object of class \code{SpatialPointsDatFrame}
#' }
#'
#' @section Package API:
#' \itemize{
#'  \item \code{\link{map_plot}} - static Base R plots
#'  \item \code{\link{map_ggplot}} - static ggplot2 plots
#'  \item \code{\link{map_ggmap}} - static ggplot2 plots with map layers
#'  \item \code{\link{map_leaflet}} - interactive Leaflet.js interactive maps
#'  \item \code{\link{map_gist}} - ineractive, shareable maps on GitHub Gists
#' }
#'
#' @importFrom methods as is
#' @importFrom stats na.omit complete.cases setNames
#' @importFrom utils write.csv browseURL
#' @importFrom graphics points
#' @importFrom ggplot2 geom_point aes ggtitle labs map_data
#' ggplot geom_point geom_polygon element_blank theme
#' @importFrom httr POST stop_for_status content upload_file
#' @importFrom sp SpatialPoints SpatialPointsDataFrame plot
#' @importFrom rworldmap getMap
#' @importFrom gistr gist_create
#' @importFrom RColorBrewer brewer.pal
#' @importFrom spocc occ2df
#' @import leaflet
#' @name mapr-package
#' @aliases mapr
#' @docType package
#' @keywords package
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
NULL
