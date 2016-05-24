pNCDun <- function(q, nu, rho, delta, n = 32, two.sided = TRUE)
{
  x <- GaussLegendre(n)
  nn <- length(q)   
  if (nn == length(nu))  xx <- cbind(q, nu) else
    if (length(nu) == 1) xx <- cbind(q, rep(nu, times = nn))  
  dtched <- function(xx)  
  {
     if (two.sided==TRUE)
     {   
        if (xx[2]=="Inf") return(pNDBD(xx[1], rho, delta, n, x)) else
           return(pNDBDF(xx[1], rho, xx[2], delta, n, x))
     } else
     {
       if (xx[2]=="Inf") return(pNDUD(xx[1], rho, delta, n, x)) else
            return(pNDUDF(xx[1], rho, xx[2], delta, n, x))
     }   
  }
  p <- apply(xx, 1, dtched)
  return(p)
}