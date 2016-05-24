additive2<-function(x,y,arg,h=1,kernel="gauss",M=2)
{
d<-length(arg)
n<-length(y)

if (kernel=="gauss") ker<-function(t){ return( exp(-t^2/2) ) }
if (kernel=="uniform") ker<-function(t){ return( (abs(t) <= 1) ) }
if (kernel=="bart") ker<-function(t){ return( (1-t) ) }

hatc<-mean(y)
residual<-matrix(y-hatc,n,1)

for (m in 1:M){
   for (j in 1:d){
      colu<-matrix(x[,j],n,1)
      estim<-matrix(0,n,d)
      if (j==1) rest<-seq(2,d) else if (j==d) rest<-seq(1:(d-1)) 
      else rest<-c(seq(1:(j-1)),seq((j+1):d))
      for (l in rest){
          for (nn in 1:n){
             curarg<-x[nn,l] 
             w<-kernesti.weights(curarg,colu,h=h,kernel=kernel)
             estim[nn,l]<-t(w)%*%residual
          }
       }
       residual<-y-hatc-matrix(rowSums(estim),n,1)
   }
}

valuevec<-matrix(0,d,1)
for (i in 1:d){
    curx<-matrix(x[,i],n,1)
    w<-kernesti.weights(arg[i],curx,h=h,kernel=kernel)
    valuevec[i]<-t(w)%*%residual
}

return(hatc+valuevec)
}


