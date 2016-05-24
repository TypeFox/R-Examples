#' ERP example data
#'
#' A data frame containing example ERP data:
#' \itemize{
#'   \item 15 subjects (S1, S2, S3, ...)
#'   \item 2 conditions (Negative and Neutral)
#' }
#' @note
#' These data have been heavily modified for use with this R package. Electrodes,
#'   subjects, and conditions have been truncated to reduce the size of the
#'   sample data set to comply with CRAN file size allowances.
#'
#' @format A data frame with 13,530 timepoints (rows) and 23
#'   columns (Subject, Stimulus, Time, electrodes)
#'
#' @source Original data collected by Hatun Zengin-Bolatkale. Data were altered
#'   for R package examples by Travis Moore.
#'
#' @examples
#' \dontrun{
#' data(ERPdata)
#' }
"ERPdata"
