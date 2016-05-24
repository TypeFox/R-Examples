
### ------ =======================================
### 12.1.1 Simulation of fractional Gaussian noise   -- book, p. 218-220
### ------ =======================================

##
## S-functions for the simulation
## of a series X(1),...,X(n) of fractional Gaussian noise
##____________________________________________________________________


ckFGN0 <- function(n,H)
{
    ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.
    if(1 > (n <- as.integer(n))) stop("'n' must be a positive integer")
    H2 <- 2*H
    k <- 0:(n-1)
    0.5 * (abs(k-1)^H2 - 2*(k)^H2 + (k+1)^H2)
}


simGauss <- function(autocov)
{
    ## Purpose: Simulation of a series X(1),...,X(n) of
    ##          a zero mean Gaussian process specified by
    ## "the first" n autocovariances  gamma(k), k=0,...,n-1

    ## ----------------------------------------------------------------------
    ## Author: Martin Maechler, Date:  5 Dec 2003, 14:48

    n <- length(autocov)
    if(n < 2) stop("need at least 2 autocovariances")

    gk <- Re(fft(autocov[1:1 + c(0:(n - 1), if(n >= 3) (n - 2):1)],
                 inverse = TRUE))

    if(any(gk < 0))
        stop("gk = Re(fft(autocov[1:n:2])) not all >= 0; invalid autocov[]")

    ## else :

    zr <- rnorm(n)
    zi <- rnorm(n-2)

    zr[c(1,n)] <- zr[c(1,n)]*sqrt(2)
    zr <- c(zr[1:n], zr[(n-1):2])
    zi <- c(0,zi,0,-zi[(n-2):1])

    z <- Re(fft((zr + 1i* zi) * sqrt(gk), inverse = TRUE))

    ts(z[1:n] / (2*sqrt(n-1)))
}

simFGN0 <- function(n,H)
{
    simGauss(ckFGN0(n,H))
}



### ------ =======================================
### 12.1.2 Simulation of fractional ARIMA($0,d,0$)   -- book, p. 220-223
### ------ =======================================

##
## Splus-functions for the simulation of a series X(1),...,X(n) of
## a fractional ARIMA(0,d,0) process  ( d = H-1/2 ).
##____________________________________________________________________

ckARMA0 <- function(n,H) {
  ## Purpose: Covariances of a fractional ARIMA(0,d,0) process
  ## -------------------------------------------------------------------------
  ## INPUT: n = length of time series
  ##        H = self-similarity parameter
  ##
  ## OUTPUT: covariances for lag  0 to n-1
  ## -------------------------------------------------------------------------
  ## Author: Jan Beran; improved: Martin Maechler, Date: Sep 95; Dec.2003

  res <- numeric(n)# k = 0:(n-1)
  d <- H - 0.5 ## d = H - 1/2  <==> 2H - 2  =  2d - 1
  g1d <- gamma(1-d)
  gd <- pi/(sin(pi*d) * g1d)# === Gamma(d): Abramowitz & Stegun 6.1.17
  res[1] <- gamma(1-2*d)/g1d^2 # == C(0) exactly
  ## Compute the first 50 ones by the exact formula:
  k <- 1:min(50,n-1)
  res[k+1] <- res[1]*gamma(k+d)* g1d/ (gamma(k-d+1)*gd)
  if(n > 51) {
      ## For large k, use the asymptotic formula -- from (2.35) :
      k <- 51:(n-1)
      res[k+1] <- res[1]* g1d/gd * k^(2*H-2)
  }
  res
}


simARMA0 <- function(n,H)
{
    ## Purpose: Simulation of a series X(1),...,X(n) of
    ##		 a fractional ARIMA(0,d,0) process (d=H-1/2)
    ## -------------------------------------------------------------------------
    ## INPUT: n = length of time series
    ##        H = self-similarity parameter
    ##
    ## OUTPUT: simulated series X(1),...,X(n)
    ## -------------------------------------------------------------------------
    ## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

    simGauss(ckARMA0(n,H))
}
## to generate  ARIMA(p,d,q), now use the above as *innovations* in
##  arima.sim(n, model=list(ar= .., ma = ..), innov= * , n.start = 0)
