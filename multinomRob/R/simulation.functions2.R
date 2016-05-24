#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1/
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#  $Id: simulation.functions2.R,v 1.3 2005/09/23 03:12:42 wrm1 Exp $
#

# Generate single random Multinomial(n,pr)
rmultinomial<-function(n=5, pr=c(0.5,0.5), long=FALSE) {
  k<-length(pr)
  if (abs(1-sum(pr))>0.000001)
   stop("(rmultinomial): parameter pr must be the k probabilities (summing to 1)")

  if(long) {
    y<-runif(n, 0, 1)
    p<-cumsum(pr)
    Seq<-1:n
    x<-sapply(y, function(y, Seq, p) {Seq[y <= p][1]}, Seq=Seq, p=p)
  } else {
    x<-rep(NA,k)
    p<-pr/c(1,(1-cumsum(pr[1:(k-1)])))
    for (i in 1:(k-1)) {
      if (n==0) {
        x[i]<-0
        if (i==k-1) x[k]<-0
        next
      }
      y<-rbinom(1,n,p[i])
      x[i]<-y
      if (i==k-1) x[k]<-n-y
      n<-n-y
    }
  }
  return(x)
}

# Generate n random Dirichlet's
rdirich<-function(n,alphavec) {
  # Gelman, et al., Bayesian Data Analysis, Chapman & Hall, 1995, p.482
  k<-length(alphavec)
  p<-matrix(rgamma(n*k,alphavec),n,k,byrow=T)
  sm<-matrix(apply(p,1,sum),n,k)
  return(p/sm)
}


