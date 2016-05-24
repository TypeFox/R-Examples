# Venkatraman E Seshan 4/20/2011
ktau <- function(x, y) {
  ok <- complete.cases(x, y)
  x1 <- x[ok]
  y1 <- y[ok]
  n <- length(x1)
  ii <- order(x1, y1)
  zzz <- .Fortran("ktau",
                  as.integer(n),
                  as.double(x1[ii]),
                  as.double(y1[ii]),
                  tau=double(1))
  zzz$tau
}
