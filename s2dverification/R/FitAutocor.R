FitAutocor <- function(estacf, window = c(-1, 1), prec = 0.01) {
  nacf <- length(estacf)
  alpha <- mean(window)
  dist <- nacf * 4
  for (jind in 1:nacf) {
    estacf[jind] <- max(estacf[jind], 0)
  }
  for (jalph in seq(window[1], window[2], prec)) {
    test <- sum(((estacf - (jalph ^ c(1:nacf - 1))) / c(1:nacf)) ^ 2)
    if (test < dist) {
      dist <- test
      alpha <- jalph
    }
  }
  alpha
}
