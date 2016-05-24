LX.mat.weibull <-
function(Y,X,sigma,phi,delta,whc) {
    r <- sum(delta)
    z <- (Y-X%*%phi)/sigma
 ndim <- ncol(X)   
 xmat <- X[,-whc]
   LX <- matrix(NA,ncol=ndim,nrow=ndim)
   
   xd <- matrix(c(xmat,z),ncol=ndim)

   for(i in 1:ndim){
       xd[,i] <- ifelse(delta==1,xd[,i],0)
         }
  LX1 <- (t(xmat)*matrix(rep(exp(z),ndim-1),nrow=nrow(t(xmat)),byrow=T))%*%xd
  LX1 <- LX1/sigma^2
  LX2 <- apply(xd/sigma^2,2,sum) - apply(matrix(rep(exp(z)*(z+1),ndim),ncol=ndim)*xd,2,sum)/sigma^2   

  LX[-ndim,] <- LX1
  LX[ndim,] <- LX2
  return(LX)
}
