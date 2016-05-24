binom.blaker.VHadj.lower.limit <- function(x,n,level,tol=1e-10,maxiter=100,nmax=n+1000,int.eps=1e-10) {
  if (x <=0) return(0)
  if (x>0) {
    lower <- binom.blaker.lower.limit(x,n,level,tol,maxiter)
    alpha <- 1-level
    nn <- n
    p <- x/n
    repeat {
      nn <- nn+1
      xstar <- nn*p
      if (qbeta(alpha/2,xstar,nn-xstar+1) > lower) break
      if (nn > nmax) {
         warning("n = ", n, ", x = ", x,": Upper limit of ", nmax, 
                 " too low, Vos-Hudson adjustment may be incomplete.") 
         break
      }
      xx <- ceiling(xstar-int.eps)
      ll <- binom.blaker.lower.limit(xx,nn,level,tol,maxiter)
      lower <- min(lower,ll)
    }
  }
  return(lower)
}
