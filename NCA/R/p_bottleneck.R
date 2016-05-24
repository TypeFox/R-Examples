p_bottleneck <-
function (loop.data, mpy, slope, intercept, cutoff, bottleneck.x) {
  bottleneck.x <- p_validate_bottleneck(bottleneck.x, "x")

  mpx <- (mpy - intercept) / slope
  if (cutoff == 0) {
    mpx [mpx < loop.data$x.low]  <- NA
    mpx [mpx > loop.data$x.high] <- NA   
  } else if (cutoff == 1) {
    mpx [mpx < loop.data$x.low]  <- loop.data$x.low
    mpx [mpx > loop.data$x.high] <- loop.data$x.high
  }

  max.value <- slope * loop.data$x.high + intercept

  # Display Xs as percentage (either cutoff or 0-high) or percentile
  if (p_bottleneck_id(bottleneck.x) == 1) {
    mpx <- 100 * (mpx - loop.data$x.low) / (loop.data$x.high - loop.data$x.low)
    max.value <- 100 * (max.value - loop.data$x.low) / (loop.data$x.high - loop.data$x.low)
  } else if (p_bottleneck_id(bottleneck.x) == 2) {
    mpx <- 100 * mpx / loop.data$x.high
    max.value <- 100 * max.value / loop.data$x.high
  } else if (p_bottleneck_id(bottleneck.x) == 4) {
    percentile <- ecdf(sort(loop.data$x))
    mpx <- matrix(100 * percentile(mpx), ncol=1)
    max.value <- 100 * percentile(max.value)
  }

  colnames(mpx) <- c(loop.data$names[1])

  return (p_pretty_mpx(loop.data, mpx, max.value))
}