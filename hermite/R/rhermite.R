rhermite <-
function(n, a, b, m=2)
{
  if (floor(m) != m ) {
    warning("improper m parameter specification")
    return(rep(NaN, n))
  }
  res <- rpois(n, a) + m*rpois(n, b)
  return(res)
}
