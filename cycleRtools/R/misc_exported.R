#' Reset a dataset or vector.
#'
#' if \code{x} is a \code{"cycleRdata"} object, all columns are reset as
#' appropriate. This can be useful after subsetting a ride dataset, for example.
#' Otherwise, this is a wrapper for \code{x - x[[1]]}.
#'
#' @param x a numeric vector or formatted cycling dataset (i.e. class \code{"cycleRdata"}).
#'
#' @return either a data frame or vector, depending on the class of \code{x}.
#'
#' @examples
#' data(ridedata)
#'
#' # Remove first minute of data and reset.
#' data_raw   <- ridedata[ridedata$timer.s > 60, ]
#' data_reset <- reset(data_raw)
#'
#' @export
reset <- function(x) UseMethod("reset", x)
#' @export
reset.default <- function(x) x - x[[1]]
#' @export
reset.cycleRdata <- function(x) {
  i    <- match(c("timer.s", "timer.min", "distance.km", "work.kJ"), colnames(x))
  x[i] <- lapply(x[i], function(y) y - y[[1]])
  x$delta.t <- c(0, Diff(x$timer.s))
  x
}
