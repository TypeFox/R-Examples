J.inf.weibul <-
function(Y,X,sigma,phi,delta,whc) {
    r <- sum(delta)
    z <- (Y-X%*%phi)/sigma
  ndim<- ncol(X)   
  xmat <- X[,-whc]
  J <- matrix(NA,ncol=ndim,nrow=ndim)
  J1 <- (t(xmat)*matrix(rep(exp(z),ndim-1),nrow=nrow(t(xmat)),byrow=T))%*%xmat
  J1 <- J1/(sigma^2) 
  J[-ndim,-ndim] <- J1
  J.LL <- sum(2*exp(z)*z/sigma^2 + (z/sigma)^2*exp(z)) - r/sigma^2 -2*sum(z/sigma^2*delta)  
  J.LP1 <- xmat*matrix(rep(exp(z),ncol(xmat)),ncol=ncol(xmat))*matrix(rep(z+1,ncol(xmat)),ncol=ncol(xmat))
  J.LP2 <- xmat*matrix(rep(delta,ncol(xmat)),ncol=ncol(xmat))
  J.LP <- (apply(J.LP1,2,sum)- apply(J.LP2,2,sum))/sigma^2
  J2 <- c(J.LP,J.LL)
  J[ndim,] <- J2
  J[,ndim] <- matrix(J2,ncol=1)
  return(J)
}
