pcf.kernesti.slice<-function(x,y,h,N,kernel="gauss",
coordi=1,p=0.5,
center=NULL,direc=NULL,radius=NULL)
{
# the slice goes through vector "vecci"
#center=rep(0,dim(x)[2]),direc=c(1,rep(0,dim(x)[2]-1)),radius=1)

n<-dim(x)[1]
d<-dim(x)[2]

if (is.null(center)){
    sl<-slice.vec(x,coordi=coordi,p=p)
    center<-sl$center
    direc<-sl$direc
    radius<-sl$radius
}

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
support[1]<--radius
support[2]<-radius
step<-(support[2]-support[1])/N
lowsuppo<-support[1]

for (i in 1:N){
     arg1d<-lowsuppo+step*i  #-step/2
     argDd<-center+arg1d*direc
     arg<-matrix(argDd,n,d,byrow=TRUE)
     w<-ker((x-arg)/h)/h^d
     w<-w/sum(w)
     value[i]<-w%*%y 
}

pcf<-list(
value=value,
down=down,high=high,
support=support,N=N)


return(pcf)
}




