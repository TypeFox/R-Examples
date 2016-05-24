#' Smorgasboard for remote sensing functions
#'
#' The package provides a variety of functions which are useful for 
#' handling, manipulating and visualizing remote sensing data.
#'
#' @name satellite-package
#' @aliases satellitepackage
#' @docType package
#' @title Smorgasboard for remote sensing functions.
#' @author Thomas Nauss, Hanna Meyer, Florian Detsch, Tim Appelhans \cr
#' \cr
#' \emph{Maintainer:} Environmental Informatics \email{admin@@environmentalinformatics-marburg.de}
#'
#' @import methods raster Rcpp
#' @importFrom stats4 plot
#' 
#' @useDynLib satellite
#' 
#' @references Some functions are taken and/or adopted from Sarah C. Goslee 
#' (2011). Analyzing Remote Sensing Data in R: The landsat Package. Journal of 
#' Statistical Software, 43(4), 1-25. Online available at 
#' \url{http://www.jstatsoft.org/v43/i04/}.
#' 
#' @keywords package
#'
NULL
#' 
#' @docType data
#' @name l7
#' @title Landsat 7 sample data
#' @description This dataset comes from the USGS. It contains part of the 
#' Landsat scene LE71950252001211EDC00 from 2001-07-30 over Maburg, Germany.
#' @details Use of this data requires your agreement to the USGS regulations on 
#' using Landsat data.
#' @format \code{raster::RasterStack} with 8 bands of 41 by 41 pixels.
#' @source
#' \url{http://earthexplorer.usgs.gov/}
NULL
#'
#' @docType data
#' @name l8
#' @title Landsat 8 sample data
#' @description This dataset comes from the USGS. It contains part of the 
#' Landsat scene LC81950252013188LGN00 from 2013-07-07 over Maburg, Germany.
#' @details Use of this data requires your agreement to the USGS regulations on 
#' using Landsat data.
#' @format \code{raster::RasterStack} with 10 bands of 41 by 41 pixels.
#' @source
#' \url{http://earthexplorer.usgs.gov/}
NULL
