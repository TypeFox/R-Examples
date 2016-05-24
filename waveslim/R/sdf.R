fdp.sdf <- function(freq, d, sigma2=1)
  sigma2 / ((2*sin(pi * freq))^2)^d

bandpass.fdp <- function(a, b, d)
  2 * integrate(fdp.sdf, lower=a, upper=b, d=d)$value

spp.sdf <- function(freq, d, fG, sigma2=1)
  sigma2 * abs(2 * (cos(2*pi*freq) - cos(2*pi*fG)))^(-2*d)

spp2.sdf <- function(freq, d1, f1, d2, f2, sigma2=1) {
  sigma2 * abs(2 * (cos(2*pi*freq) - cos(2*pi*f1)))^(-2 * d1) * 
    abs(2 * (cos(2*pi*freq) - cos(2*pi*f2)))^(-2 * d2)
}

sfd.sdf <- function(freq, s, d, sigma2=1)
  sigma2 / (2 * (1 - cos(s * 2*pi*freq)))^d

bandpass.spp <- function(a, b, d, fG) {
  if(fG > a && fG < b) {
    result1 <- integrate(spp.sdf, lower=a, upper=fG, d=d, fG=fG)$value
    result2 <- integrate(spp.sdf, lower=fG, upper=b, d=d, fG=fG)$value
  }
  else {
    result1 <- integrate(spp.sdf, lower=a, upper=b, d=d, fG=fG)$value
    result2 <- 0
  }
  return(2*(result1 + result2))
}

bandpass.spp2 <- function(a, b, d1, f1, d2, f2) {
  a1 <- a
  b1 <- b
  if(a1 < f1 && b1 > f2) {
    a2 <- f1
    b2 <- f2
    result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
    result2 <- integrate(spp2.sdf, a1, b2, d1=d1, f1=f1, d2=d2, f2=f2)$value
    result3 <- integrate(spp2.sdf, b2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
  }
  else {
    if(a1 < f1 && b1 < f2) {
      a2 <- f1
      result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
      result2 <- integrate(spp2.sdf, a2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
      result3 <- 0
    }
    else {
      if(a1 < f1 && b1 > f1 && b1 < f2) {
        a2 <- f1
        result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result2 <- integrate(spp2.sdf, a2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result3 <- 0
      }
      else {
        if(a1 > f1 && a1 < f2 && b1 > f2) {
          a2 <- f2
        result1 <- integrate(spp2.sdf, a1, a2, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result2 <- integrate(spp2.sdf, a2, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
        result3 <- 0
        }
        else {
          result1 <- integrate(spp2.sdf, a1, b1, d1=d1, f1=f1, d2=d2, f2=f2)$value
          result2 <- 0
          result3 <- 0
        }
      }
    }
  }
  return(2*(result1 + result2 + result3))
}

Hypergeometric <- function(a, b, c, z) {
  ## Recursive implementation taken from Numerical Recipes in C (6.12)
  ## Press, Teukolsky, Vetterling and Flannery (1992)
  fac <- 1
  temp <- fac
  aa <- a
  bb <- b
  cc <- c
  for(n in 1:1000) {
    fac <- fac * (aa * bb) / cc
    fac <- fac * z / n
    series <- temp + fac
    if(series == temp)
      return(series)
    temp <- series
    aa <- aa + 1
    bb <- bb + 1
    cc <- cc + 1
  }
  stop("convergence failure in Hypergeometric")
}

spp.var <- function(d, fG, sigma2=1) {
  ## Hypergeometric series representation of the variance taken from
  ## Lapsa (1997)
  omega <- 2*pi*fG
  A <- sigma2/2/sqrt(pi) * gamma(1-2*d) / gamma(3/2 - 2*d) * sin(omega)^(1-4*d)
  P1 <- Hypergeometric(1-2*d, 1-2*d, 3/2 - 2*d, sin(omega/2)^2)
  P2 <- Hypergeometric(1-2*d, 1-2*d, 3/2 - 2*d, cos(omega/2)^2)
  return(A * (P1 + P2))
}
