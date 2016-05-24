single.index<-function(x,y,arg=NULL,h=1,kernel="gauss",
M=2,method="iter",vect=FALSE,argd=arg,take=length(y),seed=1)
{
d<-dim(x)[2]

if (method=="nume") 
theta<-single.index.gene(x,y,h=h,kernel=kernel) 

if (method=="iter") 
theta<-single.index.itera(x,y,h=h,kernel=kernel,M=M,vect=vect) 

if (method=="aved") 
theta<-single.index.aved(x,y,h=h,kernel=kernel,take=take,seed=seed) 

if (method=="poid"){
   if (is.null(argd)) argd<-colMeans(x)
   theta<-single.index.aved(x,y,h=h,kernel=kernel,argd=argd) 
}

if (!is.null(arg)){
  xcur<-x%*%theta
  arg<-matrix(arg,d,1)
  acur<-sum(arg*theta)
  est<-kernesti.regr(acur,xcur,y,h=h,kernel=kernel,vect=vect)
  return(est)
}
else return(theta)
  
}




