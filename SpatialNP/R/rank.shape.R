`rank.shape` <- function(X, init=NULL, steps=Inf, eps=1e-6, maxiter=100, na.action=na.fail)
{ 
 X <- na.action(X)
 X<-as.matrix(X) 

 p<-dim(X)[2]
 if(p==1) return(diag(1))
 if (is.null(init)) init<-covshape(X)
 else init<-to.shape(init)
 if(is.finite(steps)) maxiter<-Inf

 iter<-0
 V<-init
 while(TRUE)
 {
  if(iter>=steps) return(V)
  if(iter>=maxiter) stop("maxiter reached")
  iter<-iter+1
  sqrtV<-mat.sqrt(V)
  R<-ranks(X%*%solve(sqrtV))
  V.new<-sqrtV%*%(t(R)%*%R)%*%sqrtV
  V.new<-to.shape(V.new)
  if(all(is.infinite(steps),mat.norm(V.new-V)<eps)) return(V.new)
  V<-V.new
 }
}

