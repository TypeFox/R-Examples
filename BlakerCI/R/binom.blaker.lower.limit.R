binom.blaker.lower.limit <- function(x,n,level,tol=1e-10,maxiter=100) {
  if (x <=0) return(0)
  if (x>0) {
    alpha <- 1-level
#   Clopper-Pearson limit (CPL)
    lower <- qbeta(alpha/2,x,n-x+1)
    p1 <- pbinom(x-1,n,lower,lower.tail=FALSE)
    q1.cp <- qbinom(p1,n,lower)-1
    upper <- x/n
    iter <- 0
    while (upper-lower >= tol) {
      iter <- iter+1
      if (iter > maxiter) {
        warning("Tolerance limit of ",tol, 
                             " not attained after ",maxiter, 
                             " iterations for n = ",n,", x = ",x)
        break
      }
      mid <- (lower+upper)/2
      p1 <- pbinom(x-1,n,mid,lower.tail=FALSE)
#   Blaker's limit is below the midpoint if either
#   (i)  acceptability at mid > alpha (NEW!! orig: >=), or
#   (ii) acceptability function has a discontinuity between
#        the midpoint and CPL (first test).
      if (p1 >= pbinom(q1.cp+1,n,mid) || p1 + pbinom(q1.cp,n,mid) > alpha) {
        upper <- mid
      }
      else {
        lower <- mid
      }
    }
  return(lower)
  }
}  
