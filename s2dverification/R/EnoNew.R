EnoNew <- function(xdata, detrend = FALSE, filter = FALSE) {
  alpha <- Alpha(xdata, detrend, filter)
  s <- 0
  n <- length(sort(xdata))
  for (lag in 1:n) {
    s <- s + ((n - lag) / n) * alpha ^ lag
  }

  output <- n / (1 + (2 * s))
}
