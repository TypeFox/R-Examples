simMixSNP <-
function(nSNP=5,p=0.4,ncont=1,writeFile=TRUE,outfile="sim.txt",id=1){
if(nSNP<2) stop("Error. nSNP must be  greater than 1")
if (!is.double(p) ||length(p) >1 || p < 0 || p > 1)
  stop("Allele frequency p must be a number,0<=p<=1")
if (length(ncont) >1 || ncont<1)
   stop("The number of contributors ncont must be a positive integer")
p0=p^(2*ncont)
p1=(1-p)^(2*ncont)
p2=1-p0-p1
z=sample(c(0,1,2),size=nSNP,prob=c(p0,p1,p2),replace=TRUE)
l0=length(z[z==0])
l1=length(z[z==1])
l12=length(z[z==2])

if(l0==0) z0=cbind(-1,0,0,0)
if(l0>0) z0=cbind(1:l0,rep(1,l0),rep(p,l0),rep(1,l0))

if(l1==0) z1=cbind(-1,0,0,0)
if(l1>0)z1=cbind((l0+1):(l0+l1),rep(2,l1),rep(1-p,l1),rep(1,l1))

if (l12==0) z2=cbind(-1,0,0,0)
if (l12>0) z2=cbind(sort(rep((l0+l1+1):(l0+l1+l12),2)),rep(c(1,2),l12),rep(c(p,1-p),l12),rep(1,2*l12))


simData=rbind(z0,z1,z2)
index=(1:dim(simData)[1])[simData[,1]==-1]
if(length(index)>=1) simData=simData[-index,]
d1=dim(simData)[1]
simData=data.frame(No=rep(id,d1),simData)
colnames(simData)=c("No","Marker","Allele","Frequency","Height")
if(writeFile) write.table(simData,file=outfile,quote=FALSE,row.names=FALSE)
simData
}

