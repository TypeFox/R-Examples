`penalized.pls` <-
function(X,y,P=NULL,ncomp=NULL,kernel=FALSE,scale=FALSE,blocks=1:ncol(X),select=FALSE){
n<-nrow(X)
  p<-ncol(X)

  y<-as.vector(y)

  meanx=apply(X,2,mean)

  meany=mean(y)
  if (scale==TRUE) {
    sdx = sqrt(apply(X, 2, var))
    sdx[sdx==0]=1 # take care of columns with zero variance
    }
    
  if (scale==FALSE) {
    sdx=rep(1,ncol(X))
  }
  if (is.null(ncomp)) ncomp=min(p,nrow(X)-1)
  X<-(X-rep(1,n)%*%t(meanx))/(rep(1,n)%*%t(sdx))

  y<-scale(y,center=TRUE,scale=FALSE)

   M=NULL

  if (is.null(P)==FALSE){

  Minv<-diag(p)+P

  M<-solve(Minv)
}

  if (select==FALSE){  

    if (kernel==TRUE) ppls=penalized.pls.kernel(X,y,M,ncomp); 
    if (kernel==FALSE)  ppls=penalized.pls.default(X,y,M,ncomp);
}
if (select==TRUE) ppls=penalized.pls.select(X,y,M,ncomp,blocks)
  
  
  coefficients=ppls$coefficients /(sdx %*%t(rep(1,ncol(ppls$coefficients))))
  
  
   intercept=rep(meany,ncomp) - t(coefficients)%*%meanx
  
  return(list(intercept=intercept,coefficients=coefficients))

}
