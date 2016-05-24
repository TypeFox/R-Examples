"PLACKETTsim" <-
function(n, para=NULL, ...) {
  T <- para[1]
  u <- runif(n)
  t <- runif(n)
  a <- t*(1-t)
  b <- T + a*(T-1)^2
  c <- 2*a*(u*T^2 + 1 - u) + T*(1-2*a)
  d <- sqrt(T*(T+4*a*u*(1-u)*(1-T)^2))
  v <- (c-(1-2*t)*d)/(2*b)
  z <- data.frame(U=u,V=v)
  return(z)
}
