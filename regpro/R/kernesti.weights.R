kernesti.weights<-function(arg,x,h=1,kernel="gauss",g=NULL,gernel="gauss",
vect=FALSE)
{
if (!vect) d<-length(arg) else d<-dim(x)[2]

if (d>1){

n<-dim(x)[1]

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
                      return( ans ) }
if (kernel=="exponential") ker<-function(xx){ return( exp(-rowSums(xx)) ) }

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
w<-ker((x-argu)/h)/h^d

if ((n==1)&&(is.na(w[1]))) w[1]<-1
w[is.na(w)]<-0   # make NA:s to zero

if (sum(w)==0) weights<-rep(1,n)/n else weights<-w/sum(w)

if (!is.null(g)){

   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }

   argui<-matrix(seq(n,1,-1),n,1)
   w<-ker((x-argu)/h)/h^d*ger((n-argui)/g)/g
   weights<-w/sum(w)
}
}
else{  # d==1  #########################################

n<-length(x)

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }
if (kernel=="uniform01") ker<-function(xx){ return( (abs(2*xx-1) <= 1) ) }
if (kernel=="exp") ker<-function(xx){ return( (xx>=0)*exp(-xx) ) }

if (!vect){
  x<-matrix(x,length(x),1)
  argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
  w<-ker((argu-x)/h)/h^d
  if ((n==1)&&(is.na(w[1]))) w[1]<-1
  w[is.na(w)]<-0   # make NA:s to zero
  if (sum(w)==0) weights<-rep(1,n)/n else weights<-w/sum(w)
}

if (vect){
  n<-length(x)
  x<-matrix(x,length(x),1)
  arg<-matrix(x,length(arg),1) 
  xu<-matrix(x,n,n)
  argu<-matrix(arg,n,n)
  argu<-t(argu)
  w<-ker((argu-xu)/h)/h
  w[is.na(w)]<-0   # make NA:s to zero
  #weights<-w%*%diag(1/colSums(w))
  weights<-t(t(w)/colSums(w))
}

if (!is.null(g)){

   n<-length(x)
   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }
   if (gernel=="exp") ger<-function(xx){ return( exp(-rowSums(xx))*(xx>=0) ) }

   argui<-matrix(seq(1,n,1),n,1)   # matrix(seq(n,1,-1),n,1)
   w<-ker((x-arg)/h)/h^d*ger((n-argui)/g)/g
   weights<-w/sum(w)
}

}

return(weights=weights)
}



