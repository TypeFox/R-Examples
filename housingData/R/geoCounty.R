#' County geolocation lookup table
#' 
#' Table of geographical coordinates of centroids of U.S. counties.  Computed based on a procedure described here: \url{http://stackoverflow.com/questions/9441778/improve-centering-county-names-ggplot-maps}.
#' 
#' @section Variables:
#' 
#' \itemize{
#'  \item \code{fips}: FIPS county code
#'  \item \code{county}: county name
#'  \item \code{state}: state abbreviation
#'  \item \code{lon}: county centroid longitude
#'  \item \code{lat}: county centroid latitude
#'  \item \code{rMapState}: county name as can be looked up in calls to \code{map()} in the \code{maps} package
#'  \item \code{rMapCounty}: state abbreviation as can be looked up in calls to \code{map()} in the \code{maps} package
#' }
#' @docType data
#' @name geoCounty
#' @usage geoCounty
#' @format A data frame with 3075 rows and 7 columns.
#' @examples
#' head(geoCounty)
NULL

