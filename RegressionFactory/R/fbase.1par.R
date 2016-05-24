# log-likelihood for Binomial distribution (assuming number of trials, n, is given)
# logit link function
fbase1.binomial.logit <- function(u, y, fgh=2, n=1) {
  eu <- exp(u)
  f <- -n*log(1+1/eu) - (n-y)*u # ignoring C(u) terms such as log-factorials
  if (fgh==0) return (list(f=f))
  g <- y-n/(1+1/eu)
  if (fgh==1) return (list(f=f, g=g))
  h <- -n*eu/(1+eu)^2
  return (list(f=f,g=g,h=h))
}
# probit link function
fbase1.binomial.probit <- function(u, y, fgh=2, n=1) {
  phi <- dnorm(u)
  Phi <- pnorm(u)
  f <- y*log(Phi) + (n-y)*log(1-Phi)
  if (fgh==0) return (list(f=f))
  g <- phi*(y-n*Phi)/(Phi*(1-Phi))
  if (fgh==1) return (list(f=f, g=g))
  h <- y*(u*phi/Phi - (phi/Phi)^2) - (n-y)*(u*phi/(1-Phi) + (phi/(1-Phi))^2)
  return (list(f=f,g=g,h=h))
}
# cauchit link function
fbase1.binomial.cauchit <- function(u, y, fgh=2, n=1) {
  phi <- dcauchy(u)
  Phi <- pcauchy(u)
  f <- y*log(Phi) + (n-y)*log(1-Phi)
  if (fgh==0) return (list(f=f))
  g <- y*(phi/Phi) - (n-y)*(phi/(1-Phi))
  if (fgh==1) return (list(f=f, g=g))
  h <- y*(u*phi/Phi - (phi/Phi)^2) - (n-y)*(u*phi/(1-Phi) + (phi/(1-Phi))^2)
  return (list(f=f,g=g,h=h))
}
# cloglog link function
fbase1.binomial.cloglog <- function(u, y, fgh=2, n=1) {
  eu <- exp(u)
  eeu <- exp(-eu)
  f <- y*log(1-eeu) - (n-y)*eu
  if (fgh==0) return (list(f=f))
  g <- y*(eu/(1-eeu)) - n*eu
  if (fgh==1) return (list(f=f, g=g))
  h <- eu*(y*(1-eeu-eu*eeu)/(1-eeu)^2 - n)
  return (list(f=f,g=g,h=h))
}

# log-likelihood for Poisson distribution (count response)
# log link function
fbase1.poisson.log <- function(u,y,fgh=2) {
  eu <- exp(u)
  f <- y*u-eu # ignoring additive log-factorial term since it is independent of u
  if (fgh==0) return (list(f=f))
  g <- y-eu
  if (fgh==1) return (list(f=f, g=g))
  h <- -eu
  return (list(f=f,g=g,h=h))
}

# log-likelihood for exponential distribution (positive response)
# log link function
fbase1.exponential.log <- function(u,y,fgh=2) {
  emu <- exp(-u)
  f <- -u-y*emu
  if (fgh==0) return (list(f=f))
  g <- -1+y*emu
  if (fgh==1) return (list(f=f, g=g))
  h <- -y*emu
  return (list(f=f,g=g,h=h))
}

# log-likelihood for geometric distribution (?? response)
# logit link function
# (technically, the link function applies to inverse of mean, which equals p)
fbase1.geometric.logit <- function(u,y,fgh=2) {
  eu <- exp(u)
  f <- -(y*u+(1+y)*log(1+1/eu))
  if (fgh==0) return (list(f=f))
  g <- -y+(1+y)/(1+eu)
  if (fgh==1) return (list(f=f, g=g))
  h <- -(1+y)*eu/(1+eu)^2
  return (list(f=f,g=g,h=h))
}

