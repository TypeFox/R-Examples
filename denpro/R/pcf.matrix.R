pcf.matrix<-function(A)
{
d<-2
num<-dim(A)[1]
N<-c(num,num)
recnum<-prod(N)
value<-matrix(0,recnum,1)
index<-matrix(0,recnum,d)

for (i in 1:recnum){
    inde<-digit(i-1,N)+1
    value[i]<-A[inde[1],inde[2]]
    index[i,]<-inde
}
down<-index-1
high<-index
#support<-c(0,num+1,0,num+1)
support<-c(1,num,1,num)

pcf<-list(
value=value,index=index,
down=down,high=high,  #step=delta,
support=support,N=N)
return(pcf)
}

