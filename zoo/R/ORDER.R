ORDER <- function(x, ...)
  UseMethod("ORDER")

ORDER.default <- function(x, ..., na.last = TRUE, decreasing = FALSE)
  order(x, ..., na.last = na.last, decreasing = decreasing)

ORDER.timeDate <- function(x, ...) {
  order(as.POSIXct(x), ...)
}

ORDER.chron <- ORDER.dates <- ORDER.times <- function(x, ...) {
  order(as.numeric(x), ...)
}
