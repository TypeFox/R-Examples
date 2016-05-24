#' The TM_WORLD_BORDERS_SIMPL-0.3 world map.
#' 
#' The simplified version of the world map provided by Bjorn Sandvik, thematicmapping.org.
#' 
#' The map was imported in R as follows:
#' 
#' \preformatted{
#'   require(maptools)
#'   world.map.simplified <- readShapeSpatial("~/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp")
#'   slot(world.map.simplified, 'data')[,'NAME'] <- iconv(slot(world.map.simplified, 'data')[,'NAME'], "latin1", "UTF-8")
#'   save(world.map.simplified, file="data/world.map.simplified.rda")
#' }
#' 
#' The result is a \code{SpatialPolygonsDataFrame} object. Its data slot contains a data frame with 246 observations and 11 variable:
#' 
#' \itemize{
#'   \item \strong{FIPS.} FIPS 10-4 Country Code
#'   \item \strong{ISO2.} ISO 3166-1 Alpha-2 Country Code
#'   \item \strong{ISO3.} ISO 3166-1 Alpha-3 Country Code
#'   \item \strong{UN.} ISO 3166-1 Numeric-3 Country Code
#'   \item \strong{NAME.} Name of country/area
#'   \item \strong{AREA.} Land area, FAO Statistics (2002)
#'   \item \strong{POP2005.} Population, World Polulation Prospects (2005)
#'   \item \strong{REGION.} Macro geographical (continental region), UN Statistics
#'   \item \strong{SUBREGION.} Geographical sub-region, UN Statistics
#'   \item \strong{LON.} Longitude
#'   \item \strong{LAT.} Latitude
#' }
#'
#' @note Note from the TM_WORLD_BORDERS_SIMPL-0.3's README file:
#' \itemize{
#'   \item Use this dataset with care, as several of the borders are disputed.
#'   \item The original shapefile (world_borders.zip, 3.2 MB) was downloaded from the Mapping Hacks website: http://www.mappinghacks.com/data/. The dataset was derived by Schuyler Erle from public domain sources. Sean Gilles did some clean up and made some enhancements.
#' }
#' 
#' @docType data
#' @keywords datasets
#' @format A \code{SpatialPolygonsDataFrame}.
#' @name world.map.simplified
NULL