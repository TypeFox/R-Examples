#' The time series from the M4 time series forecasting competition
#'
#' A dataset containing the prices and other attributes of almost 54,000
#' diamonds.
#'
#' @format M4 is a list (class M4data) of 10.000 objects (class M4ts) with the following structure: 
#' \describe{
#'   \item{id}{The name of the time series}
#'   \item{period}{The period of the time series. Possible values are "YEARLY", "BIANNUALLY", "QUARTERLY", "MONTHLY", "WEEKLY" & "DAILY"}
#'   \item{type}{The type of time series. Possible values are "BUSINESS-INDUSTRY", "CLIMATE", "DEMOGRAPHICS", "ECONOMICS", "FINANCE", "INTERNET-TELECOM" &"INVENTORY"}
#'   \item{H}{The forecast horizon, i.e. the number of required forecasts after the last observation}
#'   \item{n}{The number of observations in the time series}
#'   \item{past}{A time series of length n (the past observations)}
#'   \item{future}{A time series of length H (the future observations)}
#' }
#' @source \url{http://forecasters.org/resources/time-series-data/}
#'
#'@examples
#' M4
#' plot(M4[[1]])
"M4"




