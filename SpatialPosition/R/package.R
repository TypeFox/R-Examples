#' @title Spatial Position Package
#' @name SpatialPosition
#' @description Computes spatial position models: \itemize{
#' \item{Stewart potentials,}
#' \item{Reilly catchment areas,} 
#' \item{Huff catchment areas.}
#' }
#' An introduction to the package conceptual background and usage: \cr
#'  - \code{vignette(topic = "SpatialPosition")}\cr
#' A Stewart potentials use case:\cr
#'  - \code{vignette(topic = "StewartExample")}.
#' @seealso  \link{stewart}, \link{rasterStewart}, \link{plotStewart},  \link{quickStewart}, 
#' \link{rasterToContourPoly}, \link{huff}, \link{rasterHuff}, \link{plotHuff},\link{reilly}, 
#' \link{rasterReilly}, \link{plotReilly}, 
#' \link{CreateGrid}, \link{CreateDistMatrix}.
#' @docType package
NULL

#' @title Spatial Units of Paris 
#' @description A SpatialPolygonsDataFrame of the 20 spatial arrondissements of the Paris.
#' @name spatUnits
#' @docType data
NULL

#' @title Public Hospitals
#' @description A SpatialPointsDataFrame of 18 public hospitals with their capacity 
#' (Capacite field = number of beds).
#' @name spatPts
#' @docType data
NULL

#' @title Paris Perimeter
#' @description A SpatialPolygonsDataFrame of the Paris perimeter. 
#' @name spatMask
#' @docType data
NULL
