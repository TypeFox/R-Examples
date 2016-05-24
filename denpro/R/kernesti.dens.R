kernesti.dens<-function(arg,x,h=1,kernel="gauss",g=NULL,gernel="gauss")
{
d<-length(arg)

if (d>1){

if (length(h)==1) h<-rep(h,d)

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( (2*pi)^(-d/2)*exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ ans<-(rowSums(xx^2) <= 1) 
                      return( ans ) }

argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
xxx<-sweep(x-argu,2,h,"/")
w<-ker(xxx)/prod(h)
est<-sum(w)/length(w)

if (!is.null(g)){

   n<-dim(x)[1]
   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }

   argui<-matrix(seq(n,1,-1),n,1)
   w<-ker((x-argu)/h)/prod(h)*ger((n-argui)/g)/g
   est<-sum(w)/length(w)
}
}
else{  # d==1  #########################################

if (kernel=="gauss") ker<-function(xx){ return( exp(-xx^2/2) ) }
if (kernel=="uniform") ker<-function(xx){ return( (abs(xx) <= 1) ) }

x<-matrix(x,length(x),1)
w<-ker((x-arg)/h)/h^d   #weights<-w/sum(w)
est<-sum(w)/length(w)

if (!is.null(g)){

   n<-length(x)
   if (gernel=="bart") 
   ger<-function(xx){ return( (1-rowSums(xx^2))*(rowSums(xx^2)<= 1) ) }
   if (gernel=="gauss") 
   ger<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   if (gernel=="uniform") 
   ger<-function(xx){ ans<-(rowSums(xx^2)<= 1) 
                      return( ans ) }

   argui<-matrix(seq(n,1,-1),n,1)
   w<-ker((x-arg)/h)/h^d*ger((n-argui)/g)/g
   est<-sum(w)/length(w)
}

}

return(est)
}



