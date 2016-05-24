linear<-function(x,y,eleg=TRUE,lambda=0)
{
y<-matrix(y,length(y),1)
n<-dim(x)[1]
d<-dim(x)[2]

if (d==1){
  beta1<-cov(x,y)/var(x)
  beta0<-mean(y)-beta1*mean(x)
  return(list(beta0=beta0,beta1=beta1))
}
else{
  if (lambda==0){
  if (eleg){
  xnew<-matrix(0,n,d+1)
  xnew[,1]<-1
  xnew[,2:(d+1)]<-x
  alku<-t(xnew)%*%xnew
  invalku<-solve(alku,diag(rep(1,d+1)))
  beta<-invalku%*%t(xnew)%*%y
  return(list(beta0=beta[1],beta1=beta[2:(d+1)]))
  }
  else{
  kov<-cov(x)
  inkov<-solve(kov,diag(rep(1,d))) 
  beta<-inkov%*%cov(x,y)
  beta0<-mean(y)-t(beta)%*%t(t(colMeans(x)))
  return(list(beta0=beta0,beta1=beta))
  }
  }
  if (lambda>0){
  Y<-y-mean(y)
  #X<-(x-colMeans(x))/sqrt(colMeans(x^2)-colMeans(x)^2) 
  X<-x-colMeans(x)
  X<-X/sqrt(colMeans(X^2))
  XtX<-t(X)%*%X+lambda*diag(rep(1,d))
  invXtX<-solve(XtX,diag(rep(1,d)))
  beta<-invXtX%*%t(X)%*%Y
  return(list(beta0=mean(y),beta1=beta))
  }
}
}


