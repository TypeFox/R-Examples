#' SGCS: Spatial Graph based Clustering Summaries
#' 
#' Some functions for analysing clusters and clustering of spatial point patterns.
#' 
#' See package \pkg{\link{spatstat}} for more details on point pattern analysis, 
#' especially Section II: Exploratory Data Analysis.
#' 
#' This package provides functional summaries for analysing stationary and isotropic 
#' point patterns in 2D and 3D. 
#' \describe{
#' \item{ \code{\link{confun}}}{ Connectivity function, smoothed and cumulative}
#' \item{ \code{\link{clustfun}}}{ Clustering function, generalisation of clustering coefficient}
#' \item{ \code{\link{Tfun}}}{ Triangle/triplet intensity function}
#' }
#' 
#' The Ripley's \link[=Kfun]{K}-function is there as well as spatstat does not support 3D. 
#' 
#' The main source for the definitions of connectivity and clustering functions is 
#' \strong{Rajala: Spatial clustering and graph based statistical features, JYU preprints, 2010}
#' 
#' @aliases SGCS
#' @name SGCS
#' @docType package
NULL