random.pm <- function(n, values = n*diff(seq(0,1,length=n+1)**0.75)) {
  if(length(values) != n) stop("values should be of length n")
  Q <- qr.Q(qr(matrix(rnorm(n**2), nrow=n)))
  K <- Q %*% (values * t(Q))
  return( list(K = K, eigen = list(values = values, vectors = Q)) )
}
