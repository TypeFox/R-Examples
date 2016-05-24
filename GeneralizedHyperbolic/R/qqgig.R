# QQ-plot for generalised inverse gaussian distribution
qqgig <- function(y, chi = 1, psi = 1, lambda = 1,
                  param = c(chi, psi, lambda),
                  main = "GIG Q-Q Plot",
                  xlab = "Theoretical Quantiles",
                  ylab = "Sample Quantiles",
                  plot.it = TRUE, line = TRUE, ...) {

  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }

  if ((n <- length(y)) == 0)
    stop("y is empty or has only NAs")

  ## Tolerances are set very low since it doesn't matter for a plot
  x <- qgig(ppoints(n), param = param, ibfTol = 10^(-4),
            nInterpol = 201,
            uniTol = 10^(-4))[order(order(y))]

  if (has.na) {
    y <- x
    x <- yN
    x[!ina] <- y
    y <- yN
  }

  if(plot.it) {
    plot(x, y, main = main, xlab = xlab, ylab = ylab, ...)
    title(sub = paste("param = (",
          round(param[1], 3), ", " , round(param[2], 3), ", ",
          round(param[3], 3), ")", sep = ""))

    if(line)
      abline(0, 1)
  }

  invisible(list(x = x, y = y))
} ## End of qqgig()

### PP-plot for generalised inverse gaussian distribution
ppgig <- function(y, chi = 1, psi = 1, lambda = 1,
                  param = c(chi, psi, lambda),
                  main = "GIG P-P Plot",
                  xlab = "Uniform Quantiles",
                  ylab = "Probability-integral-transformed Data",
                  plot.it = TRUE, line = TRUE, ...) {

  if (has.na <- any(ina <- is.na(y))) {
    yN <- y
    y <- y[!ina]
  }

  if(0 == (n <- length(y)))
    stop("data is empty")

  ## Tolerances are set very low since it doesn't matter for a plot
  yvals <- pgig(y, param = param, ibfTol = 10^(-4))
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
          round(param[3], 3), ")", sep = ""))

    if (line)
      abline(0, 1)
  }

  invisible(list(x = xvals, y = yvals))
} ## End of ppgig()
