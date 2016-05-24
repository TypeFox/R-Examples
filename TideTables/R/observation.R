#' Sample file of high and low water times and heights
#' 
#' A sample dataset containing observation date, time and height of high and low water
#' @format A data frame with 26819 rows and 4 variables
#' \describe{
#'   \item{observation_date}{date of tide, character value in "yyyy/mm/dd" format}
#'   \item{observation_time}{time of tide, character value in "hh:mm:ss" format}
#'   \item{high_or_low_water}{indication whether high (1) or low water (0) was present given date and time., integer}
#'   \item{height}{height of the tide, numeric value}
#'}
"observation"