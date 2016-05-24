qNCDun <- function(p, nu, rho, delta, n = 32, two.sided = TRUE)
{
  x <- GaussLegendre(n)
  nn <- length(p)   
  if (nn == length(nu))  xx <- cbind(p, nu) else
    if (length(nu) == 1) xx <- cbind(p, rep(nu, times = nn))  
  dtched <- function(xx)  
  {
    if (two.sided==TRUE)
    {   
      if (xx[2]==Inf) return(qNDBD(xx[1], rho, delta, n, x)) else
        return(qNDBDF(xx[1], rho, xx[2], delta, n, x))
    } else
    {
      if (xx[2]==Inf) return(qNDUD(xx[1], rho, delta, n, x)) else
        return(qNDUDF(xx[1], rho, xx[2], delta, n, x))
    }   
  }
  d <- apply(xx, 1, dtched)
  return(d)
}