pcf.kernesti.marg<-function(x,y,h,N,kernel="gauss",coordi=1)
{
#center=rep(0,dim(x)[2]),direc=c(1,rep(0,dim(x)[2]-1)),radius=1)

n<-dim(x)[1]
d<-dim(x)[2]

if (kernel=="bart") 
   ker<-function(xx){ return( (1-rowSums(xx^2)) ) }
if (kernel=="gauss") 
   ker<-function(xx){ return( exp(-rowSums(xx^2)/2) ) }
if (kernel=="uniform") 
   ker<-function(xx){ return( (rowSums(xx^2) <= 1) ) } 

index<-seq(1:N)
len<-length(index)
down<-matrix(0,len,1)
high<-matrix(0,len,1)
down[,1]<-index-1
high[,1]<-index

value<-matrix(0,N,1)
support<-matrix(0,2,1)
support[1]<-min(x[,coordi])
support[2]<-max(x[,coordi])
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
    arg1d<-lowsuppo+step*i  #-step/2
    q<-matrix(0,1,n)
    for (ii in 1:n){
        weet<-matrix(0,n,1)
        for (j in 1:n){
            arg<-x[j,]
            arg[coordi]<-arg1d
            arg<-matrix(arg,d,1)
            w<-kernesti.weights(arg,x,h=h,kernel=kernel)
            weet[j]<-w[ii]
        }
        q[ii]<-mean(weet)
    }    
    value[i]<-q%*%y 
}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)

return(pcf)
}




