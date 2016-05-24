mstep <-
function(mu0,delta0,zm1,z0,zp1,dat,cvg,eps=1e-6){
  m <- length(mu0)
  n <- length(delta0)
  pm1 <- sum(zm1)/(m*n)
  p0 <- sum(z0)/(m*n)
  pp1 <- 1-pm1-p0
  lc1 <- function(x,y){
    a <- outer(y,x,"+")
    return(zm1*(dat-cvg*exp(a)/(1+exp(a)))+zp1*(cvg-dat-cvg*exp(a)/(1+exp(a))))
  }
  mu1 <- my.bisec(function(x) apply(lc1(x,delta0),2,sum),rep(-10,m),rep(min(-delta0),m),eps=eps)
  delta1 <- my.bisec(function(x) apply(lc1(mu1,x),1,sum),rep(-10,n),rep(min(-mu1),n),eps=eps)
  return(list(mu=mu1,delta=delta1,pRR=pm1,pRV=p0))
}
