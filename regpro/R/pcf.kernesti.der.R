pcf.kernesti.der<-function(x,y,h,N,kernel="gauss",support=NULL,direc=1,
method="ratio")
{
d<-length(N)
n<-length(y)

if (d>1){

if (kernel=="gauss"){
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
   dker<-function(xx){ 
         return( -(2*pi)^(-d/2)*xx[,direc]*exp(-rowSums(xx^2)/2) ) }
}
if (kernel=="gauss") radi<-2*h else radi<-h

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

support<-matrix(0,2*d,1)
for (i in 1:d){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
}
lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     argu<-matrix(arg,dim(x)[1],d,byrow=TRUE)
     neigh<-(rowSums((argu-x)^2) <= radi^2)
     if (sum(neigh)>=2){     # if there are obs in the neigborhood

       xred<-x[neigh,]
       yred<-y[neigh]
       argu<-matrix(arg,dim(xred)[1],d,byrow=TRUE)

       we<-ker((argu-xred)/h)/h^d
  
       w<-we/sum(we)
       u<-dker((argu-xred)/h)/h^(d+1)
       q<-1/sum(we)*(u-w*sum(u))  
  
       valli<-q%*%yred 
     }
     else valli<-mean(y)

     value[i]<-valli
     index[i,]<-inde
}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

}
else{  # d==1  #########################################

y<-matrix(y,1,length(y))
x<-matrix(x,length(x),1)

if (kernel=="gauss"){
     ker<-function(xx){ return( exp(-xx^2/2) ) }
     dker<-function(xx){ return( -xx*exp(-xx^2/2) ) }
     dker2<-function(t){ return( -t*(2*pi)^(-1/2)*exp(-t^2/2) ) }
}

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
if (is.null(support)){
   support<-matrix(0,2,1)
   support[1]<-min(x)
   support[2]<-max(x)
}
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     inde<-i
     argu<-lowsuppo+step*inde-step/2
     w<-ker((argu-x)/h)/h^1
     we<-w/sum(w)

     u<-dker((argu-x)/h)/h^(1+1)
     q<-1/sum(w)*(u-we*sum(u))  

     value[i]<-y%*%q

     if (method!="ratio"){
        xs<-sort(x)
        ys<-sort(y)
        dife<-matrix(xs[2:n]-xs[1:(n-1)],n-1,1)
        ydife<-matrix(ys[2:n],1,n-1)
        q<-dife*dker2((argu-xs[2:n])/h)/h^(1+1)
        #q<-dker2((argu-x)/h)/h^(1+1)/n
        #value[i]<-y%*%q
        value[i]<-ydife%*%q
     }
  

}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

}

return(pcf)
}
