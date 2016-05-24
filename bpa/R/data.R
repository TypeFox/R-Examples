#' Simulated Data
#'
#' Simulated (messy) data set to help illustrate some of the uses of basic
#' pattern analysis.
#'
#' \itemize{
#'   \item \code{Gender} Gender in various formats.
#'   \item \code{Date} Dates in various formats.
#'   \item{Phone} Phone numbers in various formats.
#' }
#' @docType data
#' @keywords datasets
#' @format A data frame with 1000 rows and 3 variables
#' @name messy
#' @examples
#' data(messy)
#' bpa(messy, unique_only = TRUE, ws_char = " ")
NULL
