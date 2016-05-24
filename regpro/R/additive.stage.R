additive.stage<-function(x,y=NULL,arg=NULL,residu=NULL,deet=NULL,
h=1,kernel="gauss",M=2,vect=FALSE)
{
d<-dim(x)[2]
n<-dim(x)[1]

if (is.null(residu)){
  residu<-matrix(0,n,M)
  deet<-matrix(0,M,1)
  estim<-matrix(0,n,1)
  eval<-matrix(0,n,d)
  for (m in 1:M){
    residu[,m]<-y-estim
    ssr<-matrix(0,d,1)
    estimat<-matrix(0,n,d)
    for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      ycur<-residu[,m]
      if (!vect){
         for (nn in 1:n){
             curarg<-x[nn,j]
             w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
             estimat[nn,j]<-t(w)%*%ycur
         }
      }
      else{
         estimat[,j]<-kernesti.regr(colu,colu,ycur,h=h,kernel=kernel,vect=vect)  
      }
      ssr[j]<-sum((ycur-estimat[,j])^2)
    }
    dstar<-which.min(ssr)
    deet[m]<-dstar
    eval[,dstar]<-eval[,dstar]+estimat[,dstar]
    estim<-estim+estimat[,dstar]
  }
}
else eval<-NULL

if (is.null(arg)){ 
  value<-NULL
  valuevec<-NULL
}
else{
  valuevec<-matrix(0,d,1)
  for (m in 1:M){
     ycur<-matrix(residu[,m],n,1)
     xcur<-matrix(x[,deet[m]],n,1)
     w<-kernesti.weights(arg[deet[m]],xcur,h=h,kernel=kernel)
     curvalue<-t(w)%*%ycur
     valuevec[deet[m]]<-valuevec[deet[m]]+curvalue
  }
  value<-sum(valuevec)
}

return(list(eval=eval,residu=residu,deet=deet,value=value,valvec=valuevec))
}



