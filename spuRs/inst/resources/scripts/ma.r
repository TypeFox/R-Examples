# program spuRs/resources/scripts/ma.r

ma <- function(x, w) {
  # apply a moving average of width 2*w + 1 to x
  n <- length(x)
  if (n < 2*w + 1) {
    m <- floor((n + 1)/2)
    y <- rep(0, m)
    for (i in 1:m) {
      y[i] <- mean(x[1:(2*i - 1)])
    }
  } else {
    y <- rep(0, n - w)
    for (i in 1:w) {
      y[i] <- mean(x[1:(2*i - 1)])
    }
    for (i in (w + 1):(n - w)) {
      y[i] <- mean(x[(i-w):(i+w)])
    }
  }
  return(y)
}
