#' biotic: A package for calculating a range of UK freshwater invertebrate
#' biotic indices.
#'
#' The biotic package provides a main calculation function, wrapper
#' functions for easy calculation of specific indices and a data
#' transposition function which can be used to prepare data for analysis
#' when needed.
#'
#' @section Main function:
#' The main function is \code{\link{calcindex}} which allows the
#' calculation of any of the indices implemented.
#'
#' @section Wrapper functions:
#' A function is provided for each of the individual indices to
#' allow for quick calculations. An example is
#' \code{\link{calcPSI}} which implements calculation of the
#' PSI index of sedimentation impacts.
#'
#' @section Data transposition function:
#' The \code{\link{transposedata}} function allows for simple
#' conversion between the default format with taxa in rows and
#' samples in columns and the transpose of this.
#' @importFrom stats na.omit
#' @docType package
#'
#' @name biotic-package
NULL
