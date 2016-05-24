fem.ari <- function(x,y){
  x <- as.vector(x$cls)
  y <- as.vector(y)
  xx <- outer(x, x, "==")
  yy <- outer(y, y, "==")
  upper <- row(xx) < col(xx)
  xx <- xx[upper]
  yy <- yy[upper]
  a <- sum(as.numeric(xx & yy))
  b <- sum(as.numeric(xx & !yy))
  c <- sum(as.numeric(!xx & yy))
  d <- sum(as.numeric(!xx & !yy))
  ni <- (b + a)
  nj <- (c + a)
  abcd <- a + b + c + d
  q <- (ni * nj)/abcd
  ari <- (a - q)/((ni + nj)/2 - q)
  ari
}
