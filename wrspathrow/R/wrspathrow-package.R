#' Polygon outlining TEAM site in Caxiuanã, Brazil
#' 
#' Contains a SpatialPolygonsDataFrame with a simplified polygon of the area 
#' within the Tropical Ecology Assessment and Monitoring (TEAM) network site in 
#' Caxiuanã, Brazil.
#'
#' @encoding UTF-8
#' @docType data
#' @name test_poly
NULL
.onLoad <- function(libname, pkgname) {
    load(system.file("data", "wrs1_asc_desc.RData", package="wrspathrowData"), 
         envir=parent.env(environment()))
    load(system.file("data", "wrs2_asc_desc.RData", package="wrspathrowData"),
         envir=parent.env(environment()))
}
