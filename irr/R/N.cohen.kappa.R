N.cohen.kappa<-function (rate1, rate2, k1, k0, alpha=0.05, power=0.8, twosided=FALSE) 
{
  if (twosided == FALSE) 
    d <- 1
  else d <- 2
  pi2. <- 1 - rate1
  pi.2 <- 1 - rate2
  pie <- rate1 * rate2 + pi2. * pi.2
  pi0 <- k1 * (1 - pie) + pie
  pi22 <- (pi0 - rate1 + pi.2)/2
  pi11 <- pi0 - pi22
  pi12 <- rate1 - pi11
  pi21 <- rate2 - pi11
  pi0.h <- k0 * (1 - pie) + pie
  pi22.h <- (pi0.h - rate1 + pi.2)/2
  pi11.h <- pi0.h - pi22.h
  pi12.h <- rate1 - pi11.h
  pi21.h <- rate2 - pi11.h
  Q <- (1 - pie)^(-4) * (pi11 * (1 - pie - (rate2 + rate1) * 
    (1 - pi0))^2 + pi22 * (1 - pie - (pi.2 + pi2.) * (1 - 
    pi0))^2 + (1 - pi0)^2 * (pi12 * (rate2 + pi2.)^2 + pi21 * 
    (pi.2 + rate1)^2) - (pi0 * pie - 2 * pie + pi0)^2)
  Q.h <- (1 - pie)^(-4) * (pi11.h * (1 - pie - (rate2 + rate1) * 
    (1 - pi0.h))^2 + pi22.h * (1 - pie - (pi.2 + pi2.) * 
    (1 - pi0.h))^2 + (1 - pi0.h)^2 * (pi12.h * (rate2 + pi2.)^2 + 
    pi21.h * (pi.2 + rate1)^2) - (pi0.h * pie - 2 * pie + 
    pi0.h)^2)
  N <- ((qnorm(1 - alpha/d) * sqrt(Q.h) + qnorm(power) * sqrt(Q))/(k1 - 
    k0))^2
  return(ceiling(N))
}
