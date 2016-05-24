#' Calculate The Changes In A Series
#'
#' @author George Fisher
#'
#' #@description
#'
#' @details The amount of change between two numbers in a series is current / previous.
#'
#' @export
#'
#' @examples
#' series = c(16.66, 16.85, 16.93, 16.98, 17.08, 17.03)
#' calc_chg(series)
#'
#' @return a vector of the changes with an NA in the first position
#'
#' @param series a vector of values, assumed to be an evenly-spaced time series,
#' to calculate the changes between
#'
calc_chg <- function(series) {
    chg <- series[2:length(series)] / series[1:(length(series) - 1)]
    return(c(NA, chg))
}
