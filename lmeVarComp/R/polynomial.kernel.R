polynomial.kernel <-
function(x, xi = 6L)
{
  kernel <- (1 + tcrossprod(x)) ^ xi
  kernel
}
