kernesti.der<-function(arg,x,y,h=1,direc=1,kernel="gauss",vect=FALSE)
{
d<-dim(x)[2]

if (d>1){

n<-dim(x)[1]

if (kernel=="gauss"){
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   dker<-function(xx){ 
         return( -(2*pi)^(-d/2)*xx[,direc]*exp(-rowSums(xx^2)/2) ) }
}

argu<-matrix(arg,n,d,byrow=TRUE)
we<-ker((argu-x)/h)/h^d
w<-we/sum(we)
u<-dker((argu-x)/h)/h^(d+1)
q<-1/sum(we)*(u-w*sum(u))  
value<-q%*%y 
return(value)
}

if (d==1){

if (kernel=="gauss"){
   ker<-function(xx){ return( exp(-xx^2/2) ) }
   dker<-function(xx){ return( -xx*exp(-xx^2/2) ) }
   dker2<-function(t){ return( -t*(2*pi)^(-1/2)*exp(-t^2/2) ) }
}

if (!vect){
 w<-ker((arg-x)/h)/h^1
 we<-w/sum(w)
 u<-dker((arg-x)/h)/h^(1+1)
 q<-1/sum(w)*(u-we*sum(u))  
 value<-sum(y*q)  #y%*%q
 return(value)
}
if (vect){
  n<-length(x)
  x<-matrix(x,length(x),1)
  arg<-matrix(x,length(arg),1) 
  xu<-matrix(x,n,n)
  argu<-matrix(arg,n,n)
  w<-ker((argu-xu)/h)/h^1
  we<-t(t(w)/colSums(w))
  u<-dker((argu-xu)/h)/h^(1+1)

  
  q<-1/sum(w)*(u-we*sum(u))  


  y<-matrix(y,1,n)
  value<-y%*%q
  return(value)
}

}

}



