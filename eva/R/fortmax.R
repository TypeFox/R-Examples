#' Top Ten Annual Precipitation: Fort Collins, Colorado
#'
#'
#' Top ten annual precipitation events (inches) for one rain gauge
#' in Fort Collins, Colorado from 1900 through 1999. See
#' Katz et al. (2002) Sec. 2.3.1 for more information and analyses.
#'
#' @docType data
#'
#' @usage data(fortmax)
#' @name fortmax
#'
#' @format A data frame with 100 observations. Each year is considered an observation, with the top ten annual precipitation events.
#'
#' @keywords datasets
#'
#' @references Katz, R. W., Parlange, M. B. and Naveau, P. (2002) Statistics of extremes in hydrology. Advances in Water Resources, 25, 1287-1304.
#'
#' @source Colorado Climate Center, Colorado State University. This is the original data source containing the daily precipitation data.
#'
#' @examples
#' data(fortmax)
#' y <- fortmax[, -1]
#' gevrSeqTests(y, method = "ed")
NULL
