.packageName <- 'upclass'

Aitken <- function (ll) 
{
  if (any(is.infinite(ll))) {
    linf <- Inf
    a <- NA
  }
  else {
    if (ll[2] > ll[1]) {
      a <- (ll[3] - ll[2])/(ll[2] - ll[1])
    }
    else {
      a <- 0
    }
    if (a < 1) {
      linf <- ll[3] + (ll[3] - ll[2])/(1 - a)
    }
    else {
      linf <- Inf
    }
  }
  list(ll = ll[3], linf = linf, a = a)
}
