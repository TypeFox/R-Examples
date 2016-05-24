#' @title Example event dataset
#' 
#' @description an example dataset of events, drawn from the Wikimedia web properties, for experimenting with
#' when performing session reconstruction and analysis
#'
#' @format a data.frame of 63,524 rows consisting of:
#' \describe{
#'   \item{UUID}{Hashed and salted unique identifiers representing 10,000 unique clients.}
#'   \item{timestamp}{timestamps, convertable via \code{\link{to_seconds}} (see the examples there)}
#' }
#' @name session_dataset
NULL