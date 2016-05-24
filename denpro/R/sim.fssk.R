sim.fssk<-function(n,noisedim,seed)
{
# makes n*d data matrix, d=2+noisedim
# 3 moodia, (c,0), (-c,3), (-c,-3)

d<-2+noisedim
hajo<-1
noisehajo<-sqrt(7)
c<-3^(3/2)/2
set.seed(seed)
data<-matrix(rnorm(d*n),,d)   #n*d matriisi, valkoista kohinaa
data[,1:2]<-hajo*data[,1:2]
if (noisedim>0) data[,3:d]<-noisehajo*data[,3:d]
i<-1
while (i<=n){
  mu<-matrix(0,1,d)       #moodin keskipiste
  ehto<-runif(1)
  if (ehto<1/3){          #sekoitteiden painot samat
         mu[1,1]<-0
         mu[1,2]<-c
  } 
  else if (ehto>2/3){
         mu[1,1]<-3
         mu[1,2]<--c
  }
  else{
         mu[1,1]<--3
         mu[1,2]<--c
  }
  data[i,]<-data[i,]+mu
  i<-i+1
}
return(data)
}






























