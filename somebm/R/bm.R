#' Generate a time series of Brownian motion.
#'
#' This function generatea a time series of one dimension Brownian motion,
#' adapted from http://cos.name/wp-content/uploads/2008/12/stochastic-differential-equation-with-r.pdf .
#'
#' @param x0 the start value, with the default value 0
#' @param t0 the start time point, with the default value 0
#' @param t the end time point, with the default value 1
#' @param n the number of points between t0 and t that will be generated, with the default value 100
#' @export
#' @examples
#' bm()
#' plot(bm())
#' a <- bm(x0=1, t0=1, t=2, n=1000)
#' plot(a)
bm <- function(x0=0, t0=0, t=1, n=100){
  delta <- (t-t0)/n
  W <- numeric(n+1)
  W[1] <- x0
  tseq <- seq(t0, t, length=n+1)
  for(i in 2:(n+1))
    W[i] <- W[i-1] + rnorm(1) * sqrt(delta)
  X <- ts(W, start=t0, deltat = delta)
  return(X)
}


#' Generate a time series of geometric Brownian motion.
#'
#' This function generatea a time series of one dimension geometric Brownian motion.
#' adapted from http://cos.name/wp-content/uploads/2008/12/stochastic-differential-equation-with-r.pdf .
#'
#' @param x0 the start value, with the default value 1
#' @param mu the interest rate, with the default value 0
#' @param sigma the diffusion coefficient, with the default value 1
#' @param t0 the start time point, with the default value 0
#' @param t the end time point, with the default value 1
#' @param n the number of points between t0 and t that will be generated, with the default value 100
#' @export
#' @examples
#' gbm()
#' plot(gbm())
#' b <- gbm(x0=1, mu=1, sigma=0.5, t0=1, t=2, n=1000)
#' plot(b)
gbm <- function(x0=1, mu=0, sigma=1, t0=0, t=1, n=100){
  delta <- (t-t0)/n
  W <- numeric(n+1)
  tseq <- seq(t0, t, length=n+1)
  for(i in 2:(n+1))
    W[i] <- W[i-1] + rnorm(1) * sqrt(delta)
  S <- x0 * exp((mu-sigma^2/2)*(tseq-t0) + sigma*W)
  X <- ts(S, start=t0, deltat = delta)
  return(X)
}


#' Generate a time series of fractional Brownian motion.
#'
#' This function generatea a time series of one dimension fractional Brownian motion.
#' adapted from http://www.mathworks.com.au/matlabcentral/fileexchange/38935-fractional-brownian-motion-generator .
#'
#' @param hurst the hurst index, with the default value 0.71
#' @param n the number of points between 0 and 1 that will be generated, with the default value 100
#' @export
#' @examples
#' fbm()
#' plot(fbm())
#' d <- fbm(hurst=0.2, n=1000)
#' plot(d)
fbm <- function(hurst=0.7, n=100){
  delta <- 1/n
  r <- numeric(n+1)
  r[1] <- 1
  for(k in 1:n)
    r[k+1] <- 0.5 * ((k+1)^(2*hurst) - 2*k^(2*hurst) + (k-1)^(2*hurst))
  r <- c(r, r[seq(length(r)-1, 2)])
  lambda <- Re((fft(r)) / (2*n))
  W <- fft(sqrt(lambda) * (rnorm(2*n) + rnorm(2*n)*1i))
  W <- n^(-hurst) * cumsum(Re(W[1:(n+1)]))
  X <- ts(W, start=0, deltat=delta)
  return(X)
}
