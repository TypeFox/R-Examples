S.STPPS<-function(S,x,mh)
{
S<-as.factor(S)
S<-as.factor(as.integer(S))
cum<-cumsum(mh)
sam<-matrix(0,sum(mh))
pk<-matrix(0,sum(mh))

for(k in 1: length(mh))
{
h<-which(S==k)
Nh<-length(x[h])
pkh<-x[h]/sum(x[h])
cumpk<-cumsum(pkh)
U<-runif(mh[k])
ints<-cbind(c(0,cumpk[-Nh]),cumpk)
sam.h<-rep(0,mh[k])
pk.h<-rep(0,mh[k])

for(i in 1:mh[k]){
    sam.h[i]<-which(U[i]>ints[,1] & U[i]<ints[,2])
   }
pk.h<-pkh[sam.h]

if(k==1){
sam[1:mh[k]]<-h[sam.h]
pk[1:mh[k]]<-pk.h
}
if(k>1){
sam[(cum[k-1]+1):(cum[k])]<-h[sam.h]
pk[(cum[k-1]+1):(cum[k])]<-pk.h
}

}
total<-data.frame(sam,pk)
total
}