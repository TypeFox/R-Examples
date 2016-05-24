#' Coercion Routines for Network Data Objects
#' 
#' This package contains methods for coercion between various classes used to
#' represent network data in R.
#' 
#' Functions implemented in this package allow to coerce (i.e. convert) network
#' data between classes provided by other R packages. Currently supported
#' classes are: "network" from package \pkg{network}, "igraph" from package
#' \pkg{igraph}.
#' 
#' The main functions are:
#' \itemize{
#' \item \code{\link{asNetwork}} and its methods to create objects of class
#' "network".
#' 
#' \item \code{\link{asIgraph}} and its methods to create objects of class
#' "igraph".
#' }
#' See their help pages for more information and examples.
#' 
#' As all the supported packages are written using S3 methods, so are the
#' methods in this package.
#' 
#' If you find this package useful in your work please cite it. Type
#' \code{citation(package="intergraph")} for the information how to do that.
#'
#' @importFrom network as.matrix.network
#' @importFrom network add.edges.network
#'
#' @docType package
#' @name intergraph-package
#'
#' @author Written and maintained by Michal Bojanowski \email{m.bojanowski@@icm.edu.pl}.
#'
#' @keywords package
#'
#' @example examples/package-intergraph.R
NULL
