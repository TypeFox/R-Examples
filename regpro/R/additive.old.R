additive.old<-function(x,y,arg,h=1,kernel="gauss",M=2)
{
d<-length(arg)
n<-length(y)

if (kernel=="gauss") ker<-function(t){ return( exp(-t^2/2) ) }
if (kernel=="uniform") ker<-function(t){ return( (abs(t) <= 1) ) }
if (kernel=="bart") ker<-function(t){ return( (1-t) ) }

G<-matrix(0,n,d)    # estimators g_j evaluated at x_i^j 
hatc<-mean(y)
residual<-matrix(y-hatc,n,1)

for (m in 1:M){
   for (j in 1:d){
      colu<-x[,j]
      pairdiffe<-matrix(colu,n,n,byrow=FALSE)-matrix(colu,n,n,byrow=TRUE)
      Wj<-ker(pairdiffe)
      Wj<-Wj/colSums(Wj)
      G[,j]<-t(Wj)%*%residual
      residual<-y-hatc-matrix(rowSums(G),n,1)
   }
}

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
W<-ker((x-argu)/h)/h # W<-matrix(0,n,d) kernel weights 
W<-W/colSums(W)
valuevec<-t(W)%*%residual

return(valuevec)
}




