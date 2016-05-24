# File   : pareto.r
# Author : mark van der Loo, www.markvanderloo.eu
#

# draw random sample from pareto distribution
rpareto <- function(n, xm=1, alpha=1)
{
  u <- runif(n);
  return(xm/(1-u)^(1/alpha))

}

# quantile function
qpareto <- function(p, xm=1, alpha=1)
{
  if ( alpha==0 )
   stop("alpha equals 0")
  if ( sum(p > 1 | p < 0)>0 )
   stop("p not in range [0,1]")
  
  i0 <- p==0
  q <- numeric(length(p))
  q[i0] <- xm
  q[!i0] <- xm/((1-p)^(1/alpha))
}

# density function
dpareto <- function(x, xm=1, alpha=1)
{
  I <- x > xm
  d <- 0*(1:length(x));
  d[I] <- alpha*(xm^alpha)/(x[I]^(alpha+1))
  return(d)
}

