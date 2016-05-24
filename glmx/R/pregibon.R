.ab_to_scale <- function(a, b)
  (gld::qgl(3/4, c(0, 1, a - b, a + b)) -  gld::qgl(1/4, c(0, 1, a - b, a + b)))/2.197224

.ab_to_lambda <- function(a, b)
  c(0, 1, a - b, a + b)

qpregibon <- function(p, a = 0, b = 0, lower.tail = TRUE, log.p = FALSE)
{
  if(log.p) p <- exp(p)
  if(!lower.tail) p <- 1 - p
  gld::qgl(p, .ab_to_lambda(a, b))/.ab_to_scale(a, b)
}

ppregibon <- function(q, a = 0, b = 0, lower.tail = TRUE, log.p = FALSE, tol = 1e-12)
{
  p <- gld::pgl(q * .ab_to_scale(a, b), .ab_to_lambda(a, b), inverse.eps = tol)
  if(!lower.tail) p <- 1 - p
  if(log.p) p <- log(p)
  return(p)
}

dpregibon <- function(x, a = 0, b = 0, log = FALSE, tol = 1e-12)
{
  s <- .ab_to_scale(a, b)
  d <- gld::dgl(x * s, .ab_to_lambda(a, b), inverse.eps = tol) * s
  if(log) log(d) else d
}

rpregibon <- function(n, a = 0, b = 0)
{
  qpregibon(runif(n), a = a, b = b)
}
