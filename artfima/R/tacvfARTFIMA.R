#Source: tacvfFARMA.R
#autocovariance for stationary ARTFIMA
#
"tacvfARTFIMA"   <-  
  function(d=numeric(0), lambda=numeric(0), phi = numeric(0), 
           theta = numeric(0), maxlag, sigma2 = 1)
  {
    #**WARNING**:assumes phi and theta in admissible region!!
    #use invertibleQ(phi) and invertibleQ(theta) prior to using tacvfARTFIMA!
    ARMALength <- sum(length(phi),length(theta))
    ARTFIMALength <- ARMALength+sum(length(d), length(lambda))
    isWhiteNoise <- identical(0==sum(ARTFIMALength),TRUE)
    if (isWhiteNoise) return(c(sigma2, rep(0, maxlag)))
#pure ARMA
    isARMA <- identical(length(d)==0,TRUE)||identical(d==0,TRUE) 
    if (isARMA) 
      return(tacvfARMA(phi = phi, theta = theta, maxlag = maxlag, 
                       sigma2 = sigma2))
#fractional acvf
    lagTrunc <- 2*max(128, nextn(maxlag, factors = 2))
    if (identical(lambda>1e-7,TRUE)) {
      x <- tacvfFI(d = d, lambda=lambda, maxlag = lagTrunc)
    } else {
      x <- tacvfFDWN(dfrac = d, maxlag = lagTrunc)
    }
    if (ARMALength==0) return(sigma2*x[1:(maxlag+1)])
#ARMA case
    y <- tacvfARMA(phi = phi, theta = theta, maxlag = lagTrunc, sigma2 = 1)
    z <- sigma2*mix(x, y)
    z[1:(maxlag+1)]
  }
#
#note domain of d are the real numbers
#requires package gsl
tacvfFI <-
  function(d, lambda, maxlag, sigma2=1){
    if (abs(d)<1e-8) return(c(1,rep(0,maxlag)))
    k <- 0:maxlag
    if (d>0) {
      A <- hyperg_2F1(d,d + k,1 + k,exp(-2*lambda))
      C <- k*lambda + lgamma(1+k)
      B <- lnpoch(d,k)
      ans <- A*exp(B-C)
    } else { #Approximation when d<0. When d<-0.499, acvf is neglible.
      ans <- exp(-lambda*(0:maxlag))*tacvfFDWN(max(d,-0.499), maxlag)
    }
    sigma2*ans
  }
#
"tacvfFARMA"   <-  
function(phi = numeric(0), theta = numeric(0), dfrac = 0, 
         maxlag, sigma2 = 1)
{
#phi, theta: NOT VALIDATED!!
#assumes phi and theta in admissible region and abs(d)<0.5
#white noise case
  ARMALength <- sum(length(phi),length(theta))
  isWhiteNoise <- 0==sum(ARMALength,dfrac) #works if drfrac=0 or numeric(0)
  if (isWhiteNoise) return(c(sigma2, rep(0, maxlag)))
  isARMA <- identical(0==dfrac,TRUE)||(length(dfrac)==0)
  if (isARMA) return(tacvfARMA(phi = phi, theta = theta, maxlag = maxlag, 
                          sigma2 = sigma2))
  lagTrunc <- 2*max(128, nextn(maxlag, factors = 2))
  x <- tacvfFDWN(dfrac = dfrac, maxlag = lagTrunc)
  if (ARMALength==0) return(sigma2*x[1:(maxlag+1)])
  y <- tacvfARMA(phi = phi, theta = theta, maxlag = lagTrunc, sigma2 = 1)
  z <- sigma2*mix(x, y)
  z[1:(maxlag+1)]
}
#acvf fractionally differenced white noise
#
"tacvfFDWN"   <-  
  function(dfrac, maxlag, sigma2=1)
  {
    if (dfrac>0.499) dfrac <- 0.499
    x   <-   numeric(maxlag + 1)
    x[1]   <-   gamma(1 - 2 * dfrac)/gamma(1 - dfrac)^2
    for(i in 1:maxlag)
      x[i + 1]   <-   ((i - 1 + dfrac)/(i - dfrac)) * x[i]
    x*sigma2
  }
"tacvfARMA" <-  
  function(phi = numeric(0), theta = numeric(0), maxlag = 20, 
           sigma2 = 1)
  {
    p <-length(phi)
    q <- length(theta)
    maxlagp1 <- maxlag + 1  
    if(max(p, q) == 0) {
      return(c(sigma2, numeric(maxlag)))
    }
    r   <-   max(p, q) + 1
    b   <-   numeric(r)
    C   <-   numeric(q + 1)  
    C[1]   <-   1
    theta2   <-   c(-1, theta)
    phi2   <-   numeric(3 * r)
    phi2[r]   <-   -1
    if(p > 0) {
      phi2[r + 1:p]   <-   phi
    }
    if(q > 0) {
      for(k in 1:q) {
        C[k + 1] <- (-theta[k])
        if(p > 0) {
          for(i in 1:min(p, k)) {
            C[k + 1] <- C[k + 1] + phi[i] * C[k + 1 - i]
          }
        }
      }
    }  
    for(k in 0:q) {
      for(i in k:q) {
        b[k + 1]   <-   b[k + 1] - theta2[i + 1] * C[i - k + 1]
      }
    }
    if(p == 0) {
      g   <-   c(b, numeric(maxlagp1))[1:maxlagp1]
      return(g)
    }
    else if(p > 0) {
      a   <-   matrix(numeric(r^2), ncol = r)
      for(i in 1:r) {
        for(j in 1:r) {
          if(j == 1) {
            a[i, j]   <-   phi2[r + i - 1]
          }
          else if(j != 1) {
            a[i, j] <- phi2[r + i - j] + phi2[r + i + j - 2]
          }
        }
      }
      g   <-   solve(a,  - b)
      if(length(g) <= maxlag) {
        g   <-   c(g, numeric(maxlagp1 - r)) 
        for(i in (r + 1):maxlagp1) {
          g[i]   <-   phi %*% g[i - 1:p]
        }
        return(sigma2*g[1:maxlagp1])
      }
      else if(length(g) >= maxlagp1) {
        return(sigma2*g[1:maxlagp1])
      }
    }
  }

"symtacvf"   <-  
  function(x)
  {
    c(rev(x[-1])[-1], x)
  }

mix <- function(x, y) {
  n <- 2*length(x)-2
  rev(Re(
    fft(fft(symtacvf(x)) * fft(symtacvf(y)), inverse = TRUE)/n
    )[(n/2 - 1):(n - 1)])
}




