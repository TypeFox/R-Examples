p_ceiling <- 
function (loop.data, slope, intercept) {
  if (is.na(slope) || slope < 0) {
    return(NaN)
  }

  x.low   <- loop.data$x.low
  x.high  <- loop.data$x.high
  y.low   <- loop.data$y.low
  y.high  <- loop.data$y.high

  x.ch    <- (y.high - intercept) / slope
  x.cl    <- (y.low - intercept) / slope
  y.a     <- slope * x.low  + intercept
  y.b     <- slope * x.high + intercept

  ceiling <- 0
  if (y.a > y.high) {
    if (y.b > y.high) {
      # Shouldn't be possible
      ceiling <- 0
    } else if (y.b < y.low) {
      ceiling <- (x.high - 0.5 * (x.cl+ x.ch)) * (y.high - y.low)
    } else {
      ceiling <- 0.5 * (x.high - x.ch) * (y.high - y.b)
    }
  } else if (y.a < y.low) {
    if (y.b > y.high) {
      ceiling <- (0.5 * (x.cl + x.ch) - x.low) * (y.high - y.low)
    } else if (y.b < y.low) {
      # Shouldn't be possible
      ceiling <- 0
    } else {
      ceiling <- x.low * (y.low - y.high) + x.high * (y.high - 0.5 * (y.low + y.b)) + 0.5 * x.cl * (y.b - y.low)
    }
  } else {
    if (y.b > y.high) {
      ceiling <- 0.5 * (x.ch - x.low) * (y.high - y.a)
    } else if (y.b < y.low) {
      ceiling <- 0.5 * (x.cl - x.low) * (y.a - y.low)
    } else {
      ceiling <- (x.high - x.low) * (y.high - 0.5 * (y.a + y.b))
    }
  }

  return (unname(ceiling))
}