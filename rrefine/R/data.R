#' a "dirty" data set to demonstrate rrefine features
#'
#' This data is a simulated collection of dates, days of the week, numbers of hours slept and indicators of whether or not the subject was on time for work. All observations appearing in this data set are fictitious, and any resemblance to actual arrival times for work is purely coincidental.
#'
#' @format A data frame with 63 rows and 4 variables
#' \itemize{
#'   \item {theDate} {date of observation in varying formats}
#'   \item {what.day.whas.it} {day of the week in varying formats}
#'   \item {sleephours} {number of hours slept}
#'   \item {was.i.on.time.for.work} {indicator of on-time arrival to work}
#' }
#' @examples
#' head(lateformeeting)
"lateformeeting"

#' a "clean" version of the lateformeeting sample data set
#'
#' This data is a simulated collection of dates, days of the week, numbers of hours slept and indicators of whether or not the subject was on time for work. All observations appearing in this data set are fictitious, and any resemblance to actual arrival times for work is purely coincidental.
#'
#' @format A data frame with 63 rows and 4 variables
#' \itemize{
#'   \item {date} {date of observation in POSIXct format}
#'   \item {dotw} {day of the week in consistent format}
#'   \item {hours.slept} {number of hours slept}
#'   \item {on.time} {indicator of on-time arrival to work}
#' }
#' @examples
#' head(lfm_clean)
"lfm_clean"

