##
## cramer.R
##
## Distribution of the Cramer-Von Mises test statistic
##
## $Revision: 1.2 $ $Date: 2014/06/09 04:34:49 $
##
## ..................................................................
##
##      From Matlab code written by Julian Faraway (faraway@umich.edu)
##	Translated to R by Adrian Baddeley
##
##	Reference: S. Csorgo and J.J. Faraway,
##      The exact and asymptotic distributions of Cramer-von Mises statistics
##	Journal of the Royal Statistical Society, Series B
##      58 (1996) 221-234.
##

pCvM <- local({

  ## all functions are vectorised
  D2 <- function(x) {
    z <- (x^2)/4
    b <- besselK(x=z, nu=1/4) + besselK(x=z, nu=3/4)
    b * sqrt((x^3)/(8*pi))
  }

  D3 <- function(x) {
    z <- (x^2)/4
    b <- 2*besselK(z, nu=1/4) + 3*besselK(z, nu=3/4) - besselK(z, nu=5/4)
    b * sqrt((x^5)/(32 * pi))
  }

  ED2 <- function(x) { exp(-(x^2)/4) * D2(x) }

  ED3 <- function(x) { exp(-(x^2)/4) * D3(x) }

  Ak <- function(k, x) {
    twosqrtx <- 2 * sqrt(x)
    x34 <- x^(3/4)
    x54 <- x^(5/4)
    (2*k+1)*gamma(k+1/2)*ED2((4*k+3)/twosqrtx)/(9*x34) +
      gamma(k+1/2)*ED3((4*k+1)/twosqrtx)/(72*x54) +
        2*(2*k+3)*gamma(k+3/2)*ED3((4*k+5)/twosqrtx)/(12*x54) +
          7*(2*k+1)*gamma(k+1/2)*ED2((4*k+1)/twosqrtx)/(144*x34) +
            7*(2*k+1)*gamma(k+1/2)*ED2((4*k+5)/twosqrtx)/(144*x34)
  }

  psi1 <- function(x) {
    ## Leading term in expansion of small-sample cdf of Cramer-Von Mises
    m <- length(x)
    tot <- numeric(m)
    active <- rep(TRUE, m)
    for(k in 0:20) {
      z <- -Ak(k,x[active])/(pi*factorial(k))
      tot[active] <- tot[active] + z
      active[active] <- (abs(z) >= 1e-7)
      if(!any(active)) break
    }
    return(tot + Vinf(x)/12)
  }

  Vinf <- function(x) {
    ## cdf of asymptotic distribution of Cramer-von Mises
    m <- length(x)
    tot <- numeric(m)
    active <- rep(TRUE, m)
    for(k in 0:10) {
      q <- (4*k+1)^2/(16*x[active])
      z <- ((-1)^k)*choose(-1/2,k)*sqrt(4*k+1)*
        exp(-q)*besselK(q, nu=1/4)/sqrt(x[active])
      tot[active] <- tot[active] + z
      active[active] <- (abs(z) >= 1e-7)
      if(!any(active)) break
    }
    return(tot/pi)
  }

  Vn <- function(x, n) {
    ## cdf of small-sample distribution of Cramer-von Mises statistic
    ## First order approximation, Csorgo and Faraway equation (1.8)
    Vinf(x) + psi1(x)/n
  }
    
  pCvM <- function(q, n=Inf, lower.tail=TRUE) {
    ## cdf of null distribution of Cramer-von Mises test statistic
    if(is.finite(n)) {
      lower <- 1/(12 * n)
      upper <- n/3
    } else {
      lower <- 0
      upper <- Inf
    }
    m <- length(q)
    p <- numeric(m)
    unknown <- rep(TRUE, m)
    if(any(zeroes <- (q <= lower))) {
      p[zeroes] <- 0
      unknown[zeroes] <- FALSE
    }
    if(any(ones <- (q >= upper))) {
      p[ones] <- 1
      unknown[ones] <- FALSE
    }
    if(any(unknown))
      p[unknown] <- if(is.infinite(n)) Vinf(q[unknown]) else Vn(q[unknown], n)
    return(if(lower.tail) p else 1-p)
  }

  pCvM
})

qCvM <- local({

  f <- function(x, N, P) {
    pCvM(x, N) - P
  }
    
  qCvM <- function(p, n=Inf, lower.tail=TRUE) {
    ## quantiles of null distribution of Cramer-von Mises test statistic
    stopifnot(all(p >= 0))
    stopifnot(all(p <= 1))
    if(!lower.tail) p <- 1-p
    lower <- if(is.finite(n)) (1/(12 * n)) else 0
    upper <- if(is.finite(n)) n/3 else Inf
    ans <- numeric(length(p))
    ans[p == 0] <- lower
    ans[p == 1] <- upper
    for(i in which(p > 0 & p < 1)) 
      ans[i] <- uniroot(f, c(lower, 1), N=n, P=p[i], extendInt="up")$root
    return(ans)
  }

  qCvM
})


  

