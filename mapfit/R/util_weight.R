
zero.to.inf <- list(
  phi = function(t) exp(pi*sinh(t)/2.0),
  phidash = function(t) pi*cosh(t)*exp(pi*sinh(t)/2.0)/2.0,
  interval = c(-6.8,6.8)
)

minusone.to.one <- list(
  phi = function(t) tanh(pi*sinh(t)/2.0),
  phidash = function(t) pi*cosh(t)*(1/cosh(pi*sinh(t)/2.0))^2/2.0,
  interval = c(-3,3)
)

deformula.weight <- function(f, deformula=zero.to.inf,
  zero=.Machine$double.eps, reltol=sqrt(.Machine$double.eps), ...) {
  d <- 8
  prev <- Inf
  for (k in 1:30) {
    s <- seq(from=deformula$interval[1], to=deformula$interval[2], length.out=d+1)
    h <- (deformula$interval[2] - deformula$interval[1]) / d
    res.w <- numeric(0)
    res.x <- numeric(0)
    for (t in s) {
      x <- deformula$phi(t)
      w <- deformula$phidash(t) * f(x, ...)
      if (is.finite(w) && w > zero) {
        res.w <- c(res.w, w)
        res.x <- c(res.x, x)
      }
    }
    v <- sum(res.w)*h
    aerror <- v - prev
    rerror <- aerror / prev
    if (is.finite(rerror) && abs(rerror) < reltol) {
      break
    }
    d <- d * 2
    prev <- v
  }
  list(x=res.x, weights=res.w, aerror=aerror, rerror=rerror, sum=v)
}

