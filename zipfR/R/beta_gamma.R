Cgamma <- function (a, log=!missing(base), base=exp(1))
{
  if (log) {
    if (missing(base)) {
      lgamma(a)
    } else {
      lgamma(a) / log(base)
    }
  } else {
    gamma(a)
  }
}

Rgamma <- function (a, x, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    if (missing(base)) {
      pgamma(x, shape=a, scale=1, lower.tail=lower, log.p=TRUE)
    } else {
      pgamma(x, shape=a, scale=1, lower.tail=lower, log.p=TRUE) / log(base)
    }
  } else {
    pgamma(x, shape=a, scale=1, lower.tail=lower, log.p=FALSE)
  }
}

Rgamma.inv <- function(a, y, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    if (missing(base)) {
      qgamma(y * log(base), shape=a, scale=1, lower.tail=lower, log.p=TRUE)
    } else {
      qgamma(y, shape=a, scale=1, lower.tail=lower, log.p=TRUE)
    }
  } else {
    qgamma(y, shape=a, scale=1, lower.tail=lower, log.p=FALSE)
  }
}

Igamma <- function (a, x, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    Cgamma(a, log=TRUE, base=base) + Rgamma(a, x, lower=lower, log=TRUE, base=base)
  } else {
    Cgamma(a, log=FALSE) * Rgamma(a, x, lower=lower, log=FALSE)
  }
}

Igamma.inv <- function (a, y, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    Rgamma.inv(a, y - Cgamma(a, log=TRUE, base=base), lower=lower, log=TRUE, base=base)
  } else {
    Rgamma.inv(a, y / Cgamma(a, log=FALSE), lower=lower, log=FALSE)
  }
}

Cbeta <- function(a, b, log=!missing(base), base=exp(1))
{
  if (log) {
    if (missing(base)) {
      lbeta(a, b)
    } else {
      lbeta(a, b) / log(base)
    }
  } else {
    beta(a, b)
  }
}

Rbeta <- function (x, a, b, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    if (missing(base)) {
      pbeta(x, shape1=a, shape2=b, lower.tail=lower, log.p=TRUE)
    } else {
      pbeta(x, shape1=a, shape2=b, lower.tail=lower, log.p=TRUE) / log(base)
    }
  } else {
    pbeta(x, shape1=a, shape2=b, lower.tail=lower, log.p=FALSE)
  }
}

Rbeta.inv <- function (y, a, b, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    if (missing(base)) {
      qbeta(y, shape1=a, shape2=b, lower.tail=lower, log.p=TRUE)
    } else {
      qbeta(y * log(base), shape1=a, shape2=b, lower.tail=lower, log.p=TRUE)
    }
  } else {
    qbeta(y, shape1=a, shape2=b, lower.tail=lower, log.p=FALSE)
  }
}

Ibeta <- function (x, a, b, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    Cbeta(a, b, log=TRUE, base=base) + Rbeta(x, a, b, lower=lower, log=TRUE, base=base)
  } else {
    Cbeta(a, b, log=FALSE) * Rbeta(x, a, b, lower=lower, log=FALSE)
  }
}

Ibeta.inv <- function (y, a, b, lower=TRUE, log=!missing(base), base=exp(1))
{
  if (log) {
    Rbeta.inv(y - Cbeta(a, b, log=TRUE, base=base), a, b, lower=lower, log=TRUE, base=base)
  } else {
    Rbeta.inv(y / Cbeta(a, b, log=FALSE), a, b, lower=lower, log=FALSE)
  }
}
