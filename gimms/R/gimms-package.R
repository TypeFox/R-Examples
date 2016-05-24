#' Download and Process GIMMS3g Data
#'
#' We provide a set of functions to retrieve information about GIMMS NDVI3g
#' files currently available online; download and re-arrange the bi-monthly
#' datasets according to creation time; import downloaded files from native
#' binary (ENVI) format directly into R based on the widely applied 'raster'
#' package; and calculate monthly value composites (e.g. maximum value
#' composites, MVC) from the bi-monthly input data.
#'
#' @name gimms-package
#' @aliases gimmspackage
#' @docType package
#' @title Download and Process GIMMS3g Data
#' @author Florian Detsch
#'
#' @import methods raster parallel foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom Kendall MannKendall
#' @importFrom zyp zyp.trend.vector
#' @importFrom utils download.file
#'
#' @references Pinzon, JE & Tucker, CJ (2014). A Non-Stationary 1981-2012 AVHRR
#' NDVI3g Time Series. Remote Sensing, 6(8), 6929-6960. Available online at
#' \url{http://www.mdpi.com/2072-4292/6/8/6929/htm}.
#'
#' @keywords package
#'
NULL
#'
