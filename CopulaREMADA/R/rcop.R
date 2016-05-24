# Functions for random sample bivariate copulas 
# 1-parameter copula familes with conditional cdf in closed form
# frk = Frank
# cln = Clayton (or Mardia-Takahasi-Clayton-Cook-Johnson)
# n = sample size
# cpar = copula parameter 
# Output: nx2 matrix of bivariate data with uniform(0,1) margins 

# bivariate Frank, cpar >0 or <0
rfrk=function(N,cpar)
{ u1=runif(N)
  p=runif(N)
  u2=qcondfrk(p,u1,cpar)
  cbind(u1,u2)
}


# Generates random variates from a multivariate standard normal
# distribution with linear correlation matrix x using Cholesky
# decomposition.
# Returns a matrix with n rows and siz columns. Every row is a
# random variate.
rmn<-function(n, x = matrix(c(1, 0, 0, 1), 2, 2))
{
  V <- NULL
  U <- chol(x)
  siz <- dim(x)[1]
  for(i in 1:n)
  {
    Z <- rnorm(siz)
    res <- t(U) %*% Z
    V <- cbind(V,res)
  }
  t(V)
}



rbvn=function(N,cpar)
{
  rcor<-matrix(c(1,cpar,cpar,1),2,2,byrow=T)
  dimen <- dim(rcor)[1]
  Z <- rmn(N,rcor)
  U <- NULL
  for(i in (1:N))
  U <- rbind(U, pnorm(Z[i,], 0, 1))
  U
}




# Clayton, cpar>0 
rcln=function(N,cpar)
{ u1=runif(N)
  p=runif(N)
  u2=qcondcln(p,u1,cpar)
  cbind(u1,u2)
}

rcln90=function(N,cpar)
{ cpar=-cpar
  u1=runif(N)
  p=runif(N)
  u2=qcondcln(p,u1,cpar)
  cbind(1-u1,u2)
}

rcln180=function(N,cpar)
{ u1=runif(N)
  p=runif(N)
  u2=qcondcln(p,u1,cpar)
  cbind(1-u1,1-u2)
}

rcln270=function(N,cpar)
{ cpar=-cpar
  u1=runif(N)
  p=runif(N)
  u2=qcondcln(p,u1,cpar)
  cbind(u1,1-u2)
}