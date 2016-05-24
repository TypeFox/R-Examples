#' Spatial Adjecency Data for the Polygon Grid 'shp'
#'
#' A spatial adjecency object of class 'nb' for the polygon shapefile 'shp'.
#' See documentation for \code{\link[spdep]{card}} for details on
#' structure of the object.The neighbours are according to 'queen's move' i.e
#' sharing atleast one corner. So each cell can have upto eight neibours.
#'
#' @format an object of class 'nb'
#' @name neigh_ob
NULL
