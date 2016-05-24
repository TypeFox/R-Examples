mixmat <-
function(p = 2)
{
  a <- matrix(rnorm(p * p), p, p)
  sa <- svd(a)
  d <- sort(runif(p) + 1)
  mat <- sa$u %*% (sa$v * d)
  attr(mat, "condition") <- d[p]/d[1]
  mat
}

