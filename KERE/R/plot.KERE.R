plot.KERE <- function(x,...) {
  alpha <- x$alpha[-1,]
  lambda <- x$lambda
  index <- log(lambda)
  iname <- "Log Lambda"
  xlab <- iname
  ylab <- "Coefficients"
  matplot(
        index, t(alpha), lty = 1, xlab = xlab, ylab = ylab,
        type = "l", pch = 500, 
		col = gray.colors(12, start = 0.05, end = 0.7, gamma = 2.2), ...
      )
}
