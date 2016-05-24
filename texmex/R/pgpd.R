pgpd <-
function(q, sigma, xi, u = 0, lower.tail=TRUE, log.p=FALSE ){
  ## first shift, rescale and cut at 0
  ## to handle q < u
  q <- pmax((q - u) / sigma, 0)

  ## now expand everything to the same length
  n <- max(length(q), length(xi))
  q <- rep(q, length.out=n)
  xi <- rep(xi, length.out=n)

  ## this handles negative xi properly
  xiq <- .specfun.safe.product(xi, q)

  expVal <- q * .log1prel(xiq)

  ## and now defer everything to pexp
  pexp(expVal, lower.tail=lower.tail, log.p=log.p)
}
