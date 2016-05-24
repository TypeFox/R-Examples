#' rleafmap
#'
#' @import sp
#' @import knitr
#' @import raster
#' @importFrom grDevices col2rgb dev.off heat.colors png rgb
#' @importFrom graphics par
#' @importFrom methods as is
#' @importFrom utils browseURL
#'
#' @name rleafmap
#' @docType package
NULL

#' French Hotels
#' 
#' This dataset gives the number of hotels, number of rooms and capacity for each department of metropolitan France.
#' 
#' @format a \code{SpatialPolygonsDataFrame} with geometries of the 96 french departements (epsg:4326) and 12 variables.
#' \itemize{
#'    \item DEP.CODE The code number of each department.
#'    \item DEP.NAME The name of each department.
#'    \item CHF.NAME The name of the main (administrative) city of each department.
#'    \item REGION.NAME The name of the french region (administrative) of each department.
#'    \item N.HOTELS The number of hotels.
#'    \item N.5, N.4, N.3, N.2, N.1 The number of hotels for each ranking categories (i.e. stars).
#'    \item ROOMS The number of hotel's rooms for each department.
#'    \item CAPACITY The total capacity (beds) for each department.
#'  }
#' @source \itemize{
#'    \item Institut National de l'Information Geographique et Forestiere (2014)
#'    \item ATOUT FRANCE - Agence de developpement touristique de la France (2014).
#'  }
#' @docType data
#' @keywords datasets
#' @name hotels
#' @usage data(hotels)
NULL


#' French Campsites
#' 
#' This dataset gives the number of ranked campsites and the number of tent pitches for each department of metropolitan France.
#' 
#' @format a \code{SpatialPolygonsDataFrame} with geometries of the 96 french departements (epsg:4326) and 11 variables.
#' \itemize{
#'    \item DEP.CODE The code number of each department.
#'    \item DEP.NAME The name of each department.
#'    \item CHF.NAME The name of the main (administrative) city of each department.
#'    \item REGION.NAME The name of the administrative french region of each department.
#'    \item N.CAMPSITES The number of campsites.
#'    \item N.5, N.4, N.3, N.2, N.1 The number of campsites for each ranking categories (i.e. stars).
#'    \item PITCHES The number of camp pitches for each department.
#'  }
#' @source \itemize{
#'    \item Institut National de l'Information Geographique et Forestiere (2014)
#'    \item ATOUT FRANCE - Agence de developpement touristique de la France (2014).
#'  }
#' @docType data
#' @keywords datasets
#' @name campsites
#' @usage data(campsites)
NULL


#' Velo'v stations
#' 
#' Stations of the bicycle sharing system of the city of Lyon (France).
#' 
#' @format a \code{SpatialPointsDataFrame} with the location and name of the 349 Velov stations (epsg:4326).
#' @source
#' OpenStreetMap (14/04/2014)
#'
#' @docType data
#' @keywords datasets
#' @name velov
#' @usage data(velov)
NULL
