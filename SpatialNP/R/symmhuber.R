`symmhuber` <- function(X, qg=0.9, init=NULL, steps=Inf, eps=1e-6, maxiter=100, na.action=na.fail)
{
 X <- na.action(X)
 X <- as.matrix(X) 
 p <- dim(X)[2]
 
 if (p < 2)  stop("'X' must be at least bivariate")

 if (is.null(init)) init<-cov(X)
 if(is.finite(steps)) maxiter<-Inf

 c.square <- 2 * qchisq(qg, p)
 sigma.square <- 2 * pchisq(c.square/2, p + 2) + (c.square/p) * (1 - qg)
 
 iter<-0
 V<-init

 while(TRUE)
 {
  if(iter>=steps) return(V)
  if(iter>=maxiter) stop("maxiter reached")
  iter<-iter+1
  V.new<-SSCov.hub(X,solve(V),c.square,sigma.square)
  if(all(is.infinite(steps),mat.norm(V.new-V)<eps)) 
    return(V.new)
  V<-V.new
 }
}


