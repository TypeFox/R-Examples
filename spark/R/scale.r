
#' @importFrom utils tail

scale_to <- function(data, width) {

  stopifnot(is.numeric(data), is.numeric(width),
            length(width) == 1, width == as.integer(width), width >= 2)

  sun <- seq_along(data) / length(data) * width
  sun_unit <- sun[2] - sun[1]
  res <- numeric(width)
  for (i in seq_along(res)) {
    points <- which(sun >= i-1 & sun <= i)
    w <- rep(1, length(points))

    if (points[1] != 1) {
      w[1] <- (sun[points[1]] - (i-1)) / sun_unit
    }

    if (tail(points, 1) != length(data)) {
      w1 <- (i - sun[tail(points, 1)]) / sun_unit
      w <- c(w, w1)
      points <- c(points, tail(points, 1) + 1)
    }
    res[i] <- sum(w * data[points]) / sum(w)
  }
  res
}

scale_y <- function(data, range) {
  (data - range[1]) / (range[2] - range[1])
}
