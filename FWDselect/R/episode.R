#'Episode of SO2. Pollution incident data.
#'
#'Registered values of SO2 in different
#'temporal instant. Each column of the dataset corresponds with the value
#'obtained by the series of bi-hourly means for SO2 in the instant \code{t}
#'(5-min temporal instant). The values of this dataset are greater than the
#'maximum value permitted for SO2 atmospheric.
#'
#'@name episode
#'@docType data
#'@usage episode
#'@format \code{episode} is a data frame with 19 variables (columns).
#'\describe{
#'  \item{Y}{response variable, registered values of SO2 at a specific temporal
#'  instant, in  microg/m3N. This is the value that we want to predict.}
#'  \item{In0}{registered values of SO2 at a specific temporal instant, in this
#'  case instant zero, in  microg/m3N.} \item{In1}{registered values of SO2 at a
#'  specific temporal instant, in this case 5-min instant temporal before, in
#'  microg/m3N.} \item{In2}{registered values of SO2 at a specific temporal instant,
#'  in this case 10-min instant temporal before, in microg/m3N.} ... }
#' @examples
#' data(episode)
#' head(episode)
#'


NULL
