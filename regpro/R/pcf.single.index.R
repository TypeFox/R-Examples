pcf.single.index<-function(x,y,h,N,kernel="gauss",support=NULL,
method="poid",argd=colMeans(x),type="si")
{
d<-length(N)

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

if (is.null(support)){
support<-matrix(0,2*d,1)
for (i in 1:d){
    support[2*i-1]<-min(x[,i])
    support[2*i]<-max(x[,i])
}
}

lowsuppo<-matrix(0,d,1)
for (i in 1:d) lowsuppo[i]<-support[2*i-1]
step<-matrix(0,d,1)
for (i in 1:d) step[i]<-(support[2*i]-support[2*i-1])/N[i]

if (type=="si"){

theta<-single.index(x,y,h=h,method=method,argd=argd,kernel=kernel)
xcur<-x%*%theta
for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     arg<-matrix(arg,d,1)
     acur<-sum(arg*theta)
     valli<-kernesti.regr(acur,xcur,y,h=h,kernel=kernel)
    
     value[i]<-valli
     index[i,]<-inde
}

}
else{

for (i in 1:recnum){
     inde<-digit(i-1,N)+1
     arg<-lowsuppo+step*inde-step/2
     arg<-matrix(arg,d,1)
     theta<-single.index(x,y,h=h,method="poid",argd=arg,kernel=kernel)
     acur<-sum(arg*theta) 
     xcur<-x%*%theta
     valli<-kernesti.regr(acur,xcur,y,h=h,kernel=kernel)
     value[i]<-valli
     index[i,]<-inde
}

}

down<-index-1
high<-index

pcf<-list(
value=value,index=index,
down=down,high=high,  
support=support,N=N)

return(pcf)
}


