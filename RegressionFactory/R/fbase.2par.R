# log-likelihood for Gaussian distribution, with identity and log link functions
# for mean and variance
fbase2.gaussian.identity.log <- function(u, v, y, fgh=2) {
  emv <- exp(-v)
  f <- -0.5*(v + emv*(u-y)*(u-y))
  if (fgh==0) return (list(f=f))
  g <- -cbind(emv*(u-y), 0.5*(1-emv*(u-y)*(u-y)))
  if (fgh==1) return (list(f=f, g=g))
  h <- cbind(-emv, -0.5*emv*(u-y)*(u-y), emv*(u-y))
  return (list(f=f, g=g, h=h))
}

# log-likelihood for inverse Gaussian distribution, with square
# functions for natural parameters
fbase2.inverse.gaussian.sqr.sqr <- function(u, v, y, fgh=2) {
  f <- log(abs(v)) - y*u^2 + 2*abs(u*v) - v^2/y
  if (fgh==0) return (list(f=f))
  g <- cbind(-2*y*u + 2*v, 1/v + 2*u - 2*v/y)
  if (fgh==1) return (list(f=f, g=g))
  h <- cbind(-2*y, -1/v^2-2/y, 2*sign(u)*sign(v))
  return (list(f=f, g=g, h=h))
}

fbase2.inverse.gaussian.log.log <- function(u, v, y, fgh=2) {
  emu <- exp(-u)
  emv <- exp(-v)
  f <- -0.5 * (v + emv*(y*emu - 1)^2/y)
  if (fgh==0) return (list(f=f))
  g <- cbind(emu*emv*(y*emu - 1), -0.5*(1 - emv*(y*emu - 1)^2/y))
  if (fgh==1) return (list(f=f, g=g))
  # TODO: add Hessian
}

# log-likelihood for Gamma distribution using log/log link functions for mean and dispersion parameters
# (assuming variance function of V(mu) = mu^2)
fbase2.gamma.log.log <- function(u, v, y, fgh=2) {
  emu <- exp(-u)
  emv <- exp(-v)
  logy <- log(y)
  f <- -emv*(u + y*exp(-u)) - v*emv + (emv-1)*logy - log(gamma(emv))
  if (fgh==0) return (list(f=f))
  g <- cbind(-emv*(1 - y*emu), emv*(u + y*emu - 1 + v - logy + psigamma(emv)))
  if (fgh==1) return (list(f=f, g=g))
  # TODO: add Hessian
}

