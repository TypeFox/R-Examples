direct_renorm <- function(x,y,w=NULL,sigma.sq=0) {
  n <- length(x)
  if (!is.null(w))
    W <- diag(c(sqrt(w)))
      else W <- diag(n)
  D1 <- matrix(ncol=3,nrow=n)
  D1[,1] <- x^2
  D1[,2] <- x*y
  D1[,3] <- y^2
  D1 <- W%*%D1
  D2 <- W%*%cbind(x,y,rep(1,n))
  
  S1 <- crossprod(D1)
  S2 <- crossprod(D1,D2)
  S3 <- crossprod(D2)
  
  Sx <- S3[1,3]
  Sy <- S3[2,3]
  Sxy <- S3[1,2]
  Sx2 <- S3[1,1]
  Sy2 <- S3[2,2]
  
  dS1 <- cbind(c(6*Sx2,3*Sxy,Sx2+Sy2),c(3*Sxy,Sx2+Sy2,3*Sxy),
               c(Sx2+Sy2,3*Sxy,6*Sy2))
  
  dS2 <- cbind(c(3*Sx,Sy,Sx),c(Sy,Sx,3*Sy),c(n,0,n))
  
  dS3 <- diag(c(n,n,0))
  
  S1 <- S1-sigma.sq*dS1
  S2 <- S2-sigma.sq*dS2
  S3 <- S3-sigma.sq*dS3
  
  c0<-numeric(9)
  c0[3] <- -2
  c0[5] <- 1
  c0[7]<- -2
  C1 <- matrix( c0, ncol=3)
  
  Tmatrix <- - solve(S3,t(S2))
  M <- S1 + S2 %*% Tmatrix
  M <- solve(C1,M)
  
  gev <- eigen(M)
  gevec <- gev$vectors
  a <- numeric(6)
  cond <- 4*gevec[1,]*gevec[3,]-gevec[2,]^2
  a[1:3] <- gevec[,cond > 0]
  a[4:6] <- Tmatrix %*% a[1:3]
  sig <- gev$values[cond > 0]
  
  dS <- rbind(cbind(dS1,dS2),cbind(t(dS2),dS3))
delta <- crossprod(a,dS)%*%a
  
  res <- cbind(D1,D2)%*%a
  
  return(list("residuals"=res,"fit"=list(a,gev,delta,sig)))
}
