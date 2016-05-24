additive<-function(x,y,arg=NULL,eval=NULL,
h=1,kernel="gauss",M=2,vect=FALSE)
{
d<-dim(x)[2]
n<-length(y)
hatc<-mean(y)

if (!vect){

if (is.null(eval)){
 eval<-matrix(0,n,d)
 for (m in 1:M){
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      jeval<-eval
      jeval[,j]<-0 
      ycur<-y-hatc-matrix(rowSums(jeval),n,1)
      for (nn in 1:n){
          curarg<-x[nn,j]
          w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
          eval[nn,j]<-t(w)%*%ycur
      }
      eval[,j]<-eval[,j]-mean(eval[,j])      
   }
 }
}

if (is.null(arg)){ 
  value<-NULL
  valuevec<-NULL
}
else{
  valuevec<-matrix(0,d,1)
  for (j in 1:d){
     curx<-matrix(x[,j],n,1)
     w<-kernesti.weights(arg[j],curx,h=h,kernel=kernel)
     jeval<-eval
     jeval[,j]<-0 
     ycur<-y-hatc-matrix(rowSums(jeval),n,1)
     #ycur<-matrix(eval[,j],n,1)
     valuevec[j]<-t(w)%*%ycur
  }
  value<-sum(valuevec)+hatc
}
}

######################################################
if (vect){

if (is.null(eval)){
 eval<-matrix(0,n,d)
 for (m in 1:M){
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      jeval<-eval
      jeval[,j]<-0 
      ycur<-y-hatc-matrix(rowSums(jeval),n,1)
      #############################
      xarg<-matrix(x[,j],n,1)
      W<-kernesti.weights(xarg,colu,h=h,kernel=kernel,vect=TRUE) 
      eval[,j]<-t(W)%*%ycur
      #############################
      eval[,j]<-eval[,j]-mean(eval[,j])      
   }
 }
}
if (is.null(arg)){ 
  value<-NULL
  valuevec<-NULL
}
else{
  valuevec<-matrix(0,d,1)
  for (j in 1:d){
     curx<-matrix(x[,j],n,1)
     w<-kernesti.weights(arg[j],curx,h=h,kernel=kernel)
     jeval<-eval
     jeval[,j]<-0 
     ycur<-y-hatc-matrix(rowSums(jeval),n,1)
     #ycur<-matrix(eval[,j],n,1)
     valuevec[j]<-t(w)%*%ycur
  }
  value<-sum(valuevec)+hatc
}
}

return(list(eval=eval,value=value,valvec=valuevec))
}


