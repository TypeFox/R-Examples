#' prcbench: A package to provide a testing workbench for Precision-Recall
#' curves
#'
#' The prcbench package provides four categories of important functions:
#' tool interface, test data interface, benchmarking, and curve evaluation.
#'
#' @section Tool interface:
#' The \code{\link{create_toolset}} function creates a common interface for
#'   five different tools that calculate Precision-Recall curves. These tools
#'   are \href{https://rocr.bioinf.mpi-sb.mpg.de/}{ROCR},
#'   \href{http://mark.goadrich.com/programs/AUC/}{AUCCalculator},
#'   \href{https://cran.r-project.org/package=PerfMeas}{PerfMeas},
#'   \href{https://cran.r-project.org/package=PRROC}{PRROC}, and
#'   \href{https://cran.r-project.org/package=precrec}{precrec}.
#'
#' The \code{\link{create_usrtool}} function helps users to make the same
#'  interface of the predefined ones for their own tools.

#' @section Test data interface:
#' The \code{\link{create_testset}} function creates two different types of test
#'   data sets. The first type is for benchmarking, and the second type is for
#'   curve evaluation.
#'
#' The \code{\link{create_usrdata}} function helps users to make their own test
#'   data sets.
#'
#' @section Benchmarking:
#' The \code{\link{run_benchmark}} function takes a tool set and a test data set
#'   and run \code{\link[microbenchmark]{microbenchmark}} for them.
#'
#' @section Curve evaluation:
#' The \code{\link{run_evalcurve}} function takes a tool set and a test data set
#'   and evaluates the accuracy of Precision-Recall curves for them.
#'
#' @docType package
#' @name prcbench
#'
#' @importFrom R6 R6Class
#' @importFrom ggplot2 autoplot
#' @importFrom stats runif
#' @importFrom methods slot
#' @importFrom stats aggregate
#' @importFrom methods is
#' @importFrom memoise memoise
#'
NULL

#' C1: Pre-calculated Precision-Recall curve
#'
#' A list contains scores, labels, and pre-calculated recall and precision
#' values as x and y.
#'
#' @format A list with 5 items.
#' \describe{
#'   \item{scores}{input scores}
#'   \item{labels}{input labels}
#'   \item{bp_x}{pre-calculated recall values for curve evaluation}
#'   \item{bp_y}{pre-calculated precision values for curve evaluation}
#'   \item{tp_x}{x position for displaying the test result in a plot}
#'   \item{tp_y}{y position for displaying the test result in a plot}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name C1DATA
#' @usage data(C1DATA)
NULL

#' C2: Pre-calculated Precision-Recall curve
#'
#' A list contains scores, labels, and pre-calculated recall and precision
#' values as x and y.
#'
#' @format See \code{\link{C1DATA}}.
#'
#' @docType data
#' @keywords datasets
#' @name C2DATA
#' @usage data(C2DATA)
NULL

#' C3: Pre-calculated Precision-Recall curve
#'
#' A list contains scores, labels, and pre-calculated recall and precision
#' values as x and y.
#'
#' @format See \code{\link{C1DATA}}.
#'
#' @docType data
#' @keywords datasets
#' @name C3DATA
#' @usage data(C3DATA)
NULL
