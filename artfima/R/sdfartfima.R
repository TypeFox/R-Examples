#Source: sdfArtfima.R
#spectral density function at Fourier frequencies given n
#area under sdf (-pi,pi) is 1.0
#
sdfartfima<-function(n, d, lambda, phi=numeric(0), theta=numeric(0)) {
  if (length(lambda)==0) {
    return(sdfarfima(n, d=d, phi=phi, theta=theta))
  }
  lams <- 2*pi*seq(from=1/n, to=1/2, by=1/n)
  nf <- length(lams)
  a <- outer(lams, 1:length(theta))
  if (length(theta)>0) {
    C <- cbind(1, cos(a))%*% c(1, -theta)
    S <- sin(a) %*% theta
  } else {
    C <- 1
    S <- 0
  }
  num <- as.vector(C*C + S*S)/(2*pi)#adjust area
  a <- outer(lams, 1:length(phi))
  if (length(phi)>0) {
    C <- cbind(1, cos(a))%*% c(1, -phi)
    S <- sin(a) %*% phi
  } else {
    C <- 1
    S <- 0
  }
  den <- as.vector(C*C+S*S)
  s1 <- num/den
  s2 <- (1 + exp(-2*lambda) - (2*cos(lams))/exp(lambda))^(-d)
  s1*s2
}

sdffi <- function(n, d, lambda){
  w <- 2*pi*seq(from=1/n, to=1/2, by=1/n)
  (1 + exp(-2*lambda) - (2*cos(w))/exp(lambda))^(-d)
}

sdfarfima<- function(n, d=0, phi=numeric(0), theta=numeric(0)) {
  #model assumed stationary and invertible
  lams <- 2*pi*seq(from=1/n, to=1/2, by=1/n)
  nf <- length(lams)
  a <- outer(lams, 1:length(theta))
  if (length(theta)>0) {
    C <- cbind(1, cos(a))%*% c(1, -theta)
    S <- sin(a) %*% theta
  } else {
    C <- 1
    S <- 0
  }
  num <- as.vector(C*C + S*S)/(2*pi)#adjust area
  a <- outer(lams, 1:length(phi))
  if (length(phi)>0) {
    C <- cbind(1, cos(a))%*% c(1, -phi)
    S <- sin(a) %*% phi
  } else {
    C <- 1
    S <- 0
  }
  den <- as.vector(C*C+S*S)
  s1 <- num/den
  if (!(identical(length(d)==0,TRUE) ||identical(d==0,TRUE))) {
    s2 <- (2*sin(lams/2))^(-2*d)
  }
  s1*s2
}
