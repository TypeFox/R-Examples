robacf <- function(x, lag.max = NULL, type = c("correlation", "covariance"), plot = TRUE, 
                   scaler = "s_FastQn", ...) {
  type <- match.arg(type)
  n <- length(x)

  if (is.null(lag.max))
    lag.max <- floor(10 * log10(n))
  lag.max <- min(lag.max, n - 1)

  acf <- c(1, sapply(1:lag.max, function(k) robcor(x[1:(n - k)], x[(1 + k):n], method = "ssd", scaler = scaler)))
  if (type == "covariance") {
    # acf <- acf * sapply(0:lag.max, function(k) scaler(x[1:(n - k)]) * scaler(x[(1 + k):n]))
    scaler <- match.fun(scaler)
    acf <- acf * scaler(as.numeric(x))^2
  }

  ret <- structure(
    list(
      acf = array(acf, c(lag.max + 1, 1, 1)),
      type = type,
      n.used = n,
      lag = array(c(0, 1:lag.max), c(lag.max + 1, 1, 1)),
      series = deparse(substitute(x)),
      snames = colnames(x)
    ),
    class = "acf"
  )

  if (plot) {
    plot(ret, ...)
    invisible(ret)
  } else {
    ret
  }
}
