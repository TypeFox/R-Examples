dhermite <- function (x, a, b, m = 2) 
{
  if (a < 0 || b < 0 || m < 0) {
    warning("parameters should be positive")
    return(NaN)
  }
  if (floor(m) != m || m < 2) {
    warning("improper m parameter specification")
    return(NaN)
  }
  pr <- function(x) {
    if (x%%1 != 0) {
      warning("non-integer x=", x)
      return(0)
    }
    else {
      if (a > 20 | b > 20) {
        if (edg(x, a, b, m) - edg(x - 1, a, b, m) > 0) 
          return(edg(x, a, b, m) - edg(x - 1, a, b, m))
        if (edg(x, a, b, m) - edg(x - 1, a, b, m) <= 0) 
          return(pnorm(x + 0.5) - pnorm(x - 0.5))
      }
      if (a <= 20 & b <= 20) 
        return(int.hermite(x, a, b, m))
      return(int.hermite(x, a, b, m))
    }
  }
  res <- vapply(x, pr, 1)
  return(res)
}