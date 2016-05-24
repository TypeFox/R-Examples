"dfilter" <-
function(type = c("bessel", "gauss", "custom"), param = list(pole = 4, cutoff = 1 / 10), len = ceiling(3 / param$cutoff))
{
  type <- match.arg(type)
  switch(type,
    bessel = {
      len <- as.integer(len)
      x <- 0:(len - 1)
      # filter coefficients
      param$a <- BesselPolynomial(param$pole, reverse = TRUE)
      # zeros
      param$r <- polyroot(param$a)
      # coefficients of Heaviside functions
      param$p <- sapply(1:length(param$r), function(i) 1 / prod(param$r[i] - param$r[-i]))
      # power coefficients
      A <- param$a * 1i^(0:(length(param$a) - 1))
      param$A <- sapply(1:(2 * length(A) - 1), function(i) {
        j <- max(1, i - length(A) + 1):min(i, length(A))
        sum( A[j] * Conj(A[i + 1 - j]) )
      })
      # compute cut-off frequency of "default" filter, i.e. where power is halved
      omega0 <- polyroot(param$A / param$a[1]^2 - c(2, rep(0, 2 * length(A) - 2)))
      param$omega0 <- Re(omega0[which.min(abs(Arg(omega0)))])
      # kernel function
      param$kernfun <- function(t) param$cutoff / param$omega0 * 2 * pi * Re(sapply(t * param$cutoff / param$omega0 * 2 * pi, function(s) (s > 0) * param$a[1] * sum(param$p * exp(param$r * s))))
      k <- param$kernfun(x)
      # roots and coefficients for step response
      param$rs <- c(0, param$r)
      param$ps <- sapply(1:length(param$rs), function(i) 1 / prod(param$rs[i] - param$rs[-i]))
      # step response
      param$stepfun <- function(t) Re(sapply(t * param$cutoff / param$omega0 * 2 * pi, function(s) (s > 0) * param$a[1] * sum(param$ps * exp(param$rs * s))))
      s <- param$stepfun(x)
      # power spectrum
      param$spectrum <- function(omega) Re( sapply(omega * param$omega0 / param$cutoff, function(o) param$a[1]^2 / sum(o^(1:length(param$A) - 1) * param$A) ) )
      # auto-correlation, note that the auto-correlation at lag s > 0 of Heavy(t) p exp(r t) are given by the sum over all i,j of
      # integral Heavy(t) pi exp(ri t) Heavy(t+s) pj exp(rj (t+s)) dt = integral Heavy(t) pi pj exp(rj s) exp( (ri + rj) t) = pi pj exp(rj s) / (ri + rj)
      # and for s < 0 similarly pi pj exp(ri s) / (ri + rj)
      param$acfun <- function(t) {
      acf <- sapply(c(0, t * param$cutoff / param$omega0 * 2 * pi), function(s) Re( sum( outer(1:length(param$p), 1:length(param$p), function(i, j) 
      param$p[i] * param$p[j] * exp(param$r[j] * abs(s)) / (param$r[i] + param$r[j]) )) ))
      acf[-1] / acf[1]
      }
      param$acf <- param$acfun(x)
    },
    gauss = {
      len <- as.integer(len)
      x <- 1:len - (len + 1) / 2
      k <- dnorm(x, sd = param)
      s <- pnorm(x, sd = param)
    },
    custom = {
      if(is.list(param)) {
        k <- param$kern
        s <- param$step
        param$kern <- NULL
        param$step <- NULL
      } else {
        k <- param
        s <- cumsum(k / sum(k))
      }
    }
  )
  ret <- list(type = type, param = param, kern = k / sum(k), step = s)
  # where the kernel's step response has a single jump fitted
  ret$jump <- min(which(s >= 0.5)) - 1 # last index of left half
  class(ret) <- c("dfilter", class(ret))
  ret
}

"print.dfilter" <-
function(x, ...)
{
  cat("\n")
  switch(x$type,
    bessel = {
      cat("Digitised ", x$param$pole, "-pole Bessel filter\n\n", sep = "")
      cat("cut-off frequency:", x$param$cutoff, "\n")
    },
    gauss = {
      cat("Digitised Gaussian filter\n\n")
      cat("bandwidth:", x$param, "\n")
    },
    custom = {
      cat("Custom filter\n\n", sep = "")
    }
  )
  cat("length    :", length(x$kern), "\n")
  cat("jump after:", x$jump, "\n")
  cat("step response start  :", x$step[1], "\n")
  cat("1 - step response end:", 1 - x$step[length(x$kern)], "\n")
  cat("\n")
}
