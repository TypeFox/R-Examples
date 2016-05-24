pgev <- function(q, mu, sigma, xi, lower.tail=TRUE, log.p=FALSE){
  ## first shift and scale
  q <- (q - mu) / sigma

  ## now set the lengths right
  n  <- max(length(q), length(xi))
  q  <- rep(q, length.out=n)
  xi <- rep(xi, length.out=n)

  ## this handles the limits correctly... I hope...
  xiq <- .specfun.safe.product(xi, q)

  res <- - q * .log1prel(xiq)
  pexp(exp(res), lower.tail=!lower.tail, log.p=log.p)
}


