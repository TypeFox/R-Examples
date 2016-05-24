#' @docType package
#' @name seismicRoll
#' @title Fast Rolling Statistics for Seismology
#' @description This package implements fast versions of 'roll'-ing functions primarily for use in
#' seismology. It is intended for use with the \pkg{seismic} and \pkg{seismicMetrics} packages
#' being developed for the IRIS Data Management Center (DMC) (\url{http://www.iris.edu/dms/nodes/dmc/}).
#' One advantage of the \pkg{seismicRoll} package is that all returned values are of
#' the same dimension as the incoming data with \code{NA}s where the rolling function
#' could not be calculated.
#' @details Currently exported functions include:
#' \itemize{
#'   \item{\code{\link{findOutliers}} -- outlier detection wrapper}
#'   \item{\code{\link{roll_hampel}} -- outlier detection}
#'   \item{\code{\link{roll_mean}} -- rolling mean}
#'   \item{\code{\link{roll_median}} -- rolling median (for outlier replacement)}
#'   \item{\code{\link{roll_sd}} -- rolling standard deviation}
#'   \item{\code{\link{roll_stalta}} -- first break picker}
#' }
#' 
#' \strong{History}
#' 
#' version 1.0.0 -- initial release
#' 
#' @useDynLib seismicRoll
#' @importFrom Rcpp sourceCpp
NULL
