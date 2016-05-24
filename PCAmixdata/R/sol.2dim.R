sol.2dim <- function(A,indexj,p,p1)
{
  uj <- rep(NA,p)
  vj <- rep(NA,p)
  v1 <- A[,1]^2-A[,2]^2
  v2 <- A[,1]*A[,2]
  if (p1 >0) {
    uj[1:p1]<-v1[1:p1]
    vj[1:p1]<-2*v2[1:p1]
  }
  if ((p-p1) >0) 
    for (j in (p1+1):p) {
      uj[j] <- sum(v1[which(indexj==j)])
      vj[j] <- 2*sum(v2[which(indexj==j)])
    }
  a <- 2*p*sum(uj*vj) - 2*sum(uj)*sum(vj)
  b <- p*sum(uj^2-vj^2) - sum(uj)^2+sum(vj)^2
  
  if (a>=0) {
    r <- acos(b/sqrt(a^2+b^2))
    theta <- r/4
  } else  {
    r <- -acos(b/sqrt(a^2+b^2))
    theta <- r/4
  }
  
  theta <- as.numeric(theta)
  T <- matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),2,2,byrow=TRUE)
  return(list(theta=theta,T=T))
  
}