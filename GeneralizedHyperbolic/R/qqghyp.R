### Q-Q plot for generalized hyperbolic distribution
qqghyp <- function(y, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                   param = c(mu, delta, alpha, beta, lambda),
                   main = "Generalized Hyperbolic Q-Q Plot",
                   xlab = "Theoretical Quantiles",
                   ylab = "Sample Quantiles",
                   plot.it = TRUE, line = TRUE, ...) {

  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }

  if ((n <- length(y)) == 0)
    stop("y is empty or has only NAs")

  x <- qghyp(ppoints(n), param = param)[order(order(y))]

  if (has.na) {
    y <- x
    x <- yN
    x[!ina] <- y
    y <- yN
  }

  if (plot.it) {
    plot(x, y, main = main, xlab = xlab, ylab = ylab, ...)
    title(sub = paste("param = (",
            round(param[1], 3), ", ", round(param[2], 3), ", ",
            round(param[3], 3), ", ", round(param[4], 3), ", ",
            round(param[5], 3), ")", sep = ""))

    if (line)
      abline(0, 1)
  }

  invisible(list(x = x, y = y))
} ## End of qqghyp()

### P-P plot for generalized hyperbolic distribution
ppghyp <- function(y, mu = 0, delta = 1, alpha = 1, beta = 0, lambda = 1,
                   param = c(mu, delta, alpha, beta, lambda),
                   main = "Generalized Hyperbolic P-P Plot",
                   xlab = "Uniform Quantiles",
                   ylab = "Probability-integral-transformed Data",
                   plot.it = TRUE, line = TRUE, ...) {

  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }

  if((n <- length(y)) == 0)
    stop("data is empty")

  yvals <- pghyp(y, param = param)
  xvals <- ppoints(n, a = 1 / 2)[order(order(y))]

  if (has.na) {
    y <- yvals
    x <- xvals
    yvals <- yN
    yvals[!ina] <- y
    xvals <- yN
    xvals[!ina] <- x
  }

  if (plot.it) {
    plot(xvals, yvals, main = main, xlab = xlab, ylab = ylab,
         ylim = c(0, 1), xlim = c(0, 1), ...)
    title(sub = paste("param = (",
            round(param[1], 3), ", ", round(param[2], 3), ", ",
            round(param[3], 3), ", ", round(param[4], 3), ", ",
            round(param[5], 3), ")", sep = ""))

    if (line)
      abline(0, 1)
  }

  invisible(list(x = xvals, y = yvals))
} ## End of ppghyp()
