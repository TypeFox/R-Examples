#' Emission of SO2. Pollution incident data.
#'
#' Registered values of SO2 in different temporal instant. Each column of the
#' dataset corresponds with the value obtained by the series of bi-hourly means
#' for SO2 in the instant \eqn{t} (5-min temporal instant).
#'
#'@name pollution
#' @docType data
#' @usage pollution
#' @format \code{pollution} is a data frame with 19 variables (columns).
#'   \describe{
#'   \item{Y}{response variable, registered values of SO2 at a
#'   specific temporal instant, in  microg/m3N. This is the value that we want to
#'   predict.}
#'   \item{In0}{registered values of SO2 at a specific temporal
#'   instant, in this case instant zero, in  microg/m3N.}
#'   \item{In1}{registered
#'   values of SO2 at a specific temporal instant, in this case 5-min instant
#'   temporal before, in  microg/m3N.}
#'  \item{In2}{registered values of SO2 at a
#'   specific temporal instant, in this case 10-min instant temporal before, in
#'   microg/m3N.}
#'   ...
#'   }
#' @examples
#' data(pollution)
#' head(pollution)
#'


NULL
