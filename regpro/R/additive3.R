additive3<-function(x,y,arg,h=1,kernel="gauss",M=2)
{
d<-length(arg)
n<-length(y)

hatc<-mean(y)
estim<-matrix(0,n,d)
for (m in 1:M){
   itestim<-estim
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      for (nn in 1:n){
          curarg<-x[nn,j]
          jestim<-itestim
          jestim[,j]<-0 
          ycur<-y-hatc-matrix(rowSums(jestim),n,1)
          w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
          itestim[nn,j]<-t(w)%*%ycur
      }
      itestim[,j]<-itestim[,j]-mean(itestim[,j])      
   }
   estim<-itestim
}

valuevec<-matrix(0,d,1)
for (j in 1:d){
    curx<-matrix(x[,j],n,1)
    w<-kernesti.weights(arg[j],curx,h=h,kernel=kernel)
    jestim<-estim
    jestim[,j]<-0 
    ycur<-y-hatc-matrix(rowSums(jestim),n,1)
    valuevec[j]<-t(w)%*%ycur
}

return(estim)
#return(hatc+valuevec)
}



