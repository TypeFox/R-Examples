.yw <- function(r, k = length(r) - 1, sigma.too = FALSE) {
  mat <- toeplitz(r[1:k])
  stopifnot(all(eigen(mat, T, T)$values > 0))
  phi <- solve(mat, r[2:(k + 1)])
  if (sigma.too) {
    sigma <- r[1] - sum(phi * r[2:(k + 1)])
    list(phi = phi, sigma = sigma)
  } else {
    phi
  }
}

robar <- function(x, order = 2, scaler = "s_FastQn") {
  r <- drop(robacf(x, lag.max = order, type = "cov", plot = FALSE, scaler = scaler)$acf)
  ar <- .yw(r, k = order, sigma.too = TRUE)

  ret <- structure(
    list(
      order=order,
      ar=ar$phi,
      var.pred=ar$sigma,
      n.used=length(x),
      method = "Robust Yule-Walker",
      series=deparse(substitute(x)),
      frequency=frequency(x),
      call=match.call()
    ),
    class = "ar"
  )
  ret
}

# example:
#
# x <- arima.sim(list(ar=c(1,-0.7)), 100)
# y <- drop(robacf(x)$acf)
# phi <- .yw(y, 2)
# s <- function(phi, f) 1/(1 + phi[1]^2 + phi[2]^2 - 2 * phi[1] * (1 - phi[2]) * cos(2*pi*f) - 2 * phi[2] * cos(4*pi*f))
# sx <- (0:500)/1000
# sy <- sapply(sx, function(f) s(phi, f))
# plot(sx, sy)
#
# spec.ar(robar(x, 2))
# spec.ar(robar(ldeaths, 10))
