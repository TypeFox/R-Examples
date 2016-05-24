################################
#### Hypothesis test for a mean direction
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, pg. 212
################################

meandir.test <- function(x, mu, B = 999) {
  ## x is the sample
  ## mu is the hypothesized mean direction under H0
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x/sqrt( rowSums(x^2) )  ## makes sure x are unit vectors
  mu <- mu/sqrt( sum(mu^2) )  ## makes sure m0 is a unit vector
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  k1 <- vmf(x)$k  ## concentration parameter under H1
  xbar <- colMeans(x)  ## x-bar
  m1 <- xbar/sqrt(sum(xbar^2))
  lik <- function(k) {
    n * (p/2 - 1) * log(k) - 0.5 * n * p * log(2 * pi) + k * sum(x %*% mu) -
      n * (log(besselI(k, p/2 - 1, expon.scaled = T)) + k)
  }
  qa0 <- optimize(lik, c(0, 100000), maximum = T)  ## log-likelihood under H0
  k0 <- qa0$maximum  ## concentration parameter under H0
  apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) +
  log(besselI(k0, p/2 - 1, expon.scaled = T)) + k0
  apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) +
  log(besselI(k1, p/2 - 1, expon.scaled = T)) + k1
  w <- 2 * n * (k1 * sqrt(sum(xbar^2)) - k0 * sum(mu * xbar) - apk1 + apk0)
  if (B == 1) {
    pvalue <- 1 - pchisq(w, p - 1)
  }
  if (B > 1) {
    A <- rotation(m1, mu)
    y <- x %*% t(A)  ## bring the data under H0
    ## y has mean direction equal to mu
    wb <- numeric(B)
    for (i in 1:B) {
      nu <- sample(1:n, n, replace = T)
      z <- y[nu, ]
      k1 <- vmf(z)$k  ## concentration parameter under H1
      zbar <- colMeans(z)  ## z-bar
      lik <- function(k) {
       n * (p/2 - 1) * log(k) - 0.5 * n * p * log(2 * pi) +
       k * sum(z %*% mu) - n * (log(besselI(k, p/2 - 1, expon.scaled = TRUE)) + k)
      }
      qa0 <- optimize(lik, c(0, 100000), maximum = TRUE)  ## log-likelihood under H0
      k0 <- qa0$maximum  ## concentration parameter under H0
      apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) +
      log(besselI(k0, p/2 - 1, expon.scaled = TRUE)) + k0
      apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) +
      log(besselI(k1, p/2 - 1, expon.scaled = TRUE)) + k1
      wb[i] <- 2 * n * (k1 * sqrt(sum(zbar^2)) - k0 * sum(mu * zbar) -
      apk1 + apk0)
    }
    pvalue <- (sum(wb > w) + 1)/(B + 1)
  }
  list(mean.dir = m1, pvalue = pvalue)
}
