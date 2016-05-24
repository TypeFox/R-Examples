#' The \pkg{ezec} package for easy EC calculation.
#'
#' @section Introduction: The package \pkg{ezec} is not a revolutionary work. It
#'   simply is a wrapper for the \pkg{drc} package that makes life a little
#'   easier when it comes to calculating a simple EC 50. The main function of
#'   the package is \code{\link{EC_table}}. This function will do as it says and
#'   automatically produce a table to EC values for each isolate in your sample.
#'
#' @section Data format: Data is expected to exist in a table with at least
#'   three columns: \itemize{ \item Sample ID \item Dosage \item Response value
#'   (Growth) } Any other columns in your data are optional. An example data set
#'   is \code{\link{dummydata}}.
#'
#' @name ezec
#' @docType package
#' @import dplyr drc
#' @importFrom graphics plot plot.new text
#' @importFrom stats na.omit
#' @importFrom utils read.table
NULL

#' dummydata
#'
#' @name dummydata
#' @docType data
#' @usage data(dummydata)
#' @format a data frame with 96 rows and 7 columns representing two isolates
#'   tested for Metalaxyl resistance over 6 concentrations with 8 replicates per
#'   concentration. Each rep number were conducted in separate weeks. The First
#'   sample is real and the second is fake.
NULL
