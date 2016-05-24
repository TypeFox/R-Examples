Qfun <- function(theta, suff.stat, n) {
  mu<-rep(0,2)
  Sigma<-matrix(0, 2,2)
  Suff1<-rep(0,2)
  Suff2<-matrix(0,2,2)
  
  mu <- theta[1:2]
  Sigma[1,1]<-theta[3]
  Sigma[2,2]<-theta[4]
  Sigma[1,2]<-Sigma[2,1]<-theta[5]*sqrt(Sigma[1,1]*Sigma[2,2])

  Suff1 <- n*suff.stat[1:2]
  Suff2[1,1]<-n*suff.stat[3]
  Suff2[2,2]<-n*suff.stat[4]
  Suff2[1,2]<-n*suff.stat[5]
 
  invSigma<-solve(Sigma)

  return(-0.5*n*log(det(Sigma))-0.5*sum(diag(invSigma%*%(Suff2-mu%*%t(Suff1)-Suff1%*%t(mu)+n*mu%*%t(mu)))))

}
