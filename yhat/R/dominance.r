dominance <- function (apsOut){

varid<-apsOut$ivID
PredBitMap<-apsOut$PredBitMap
indin<-ind<-apsOut$apsBitMap
ind[indin[,1],1]<-1:nrow(indin)
APSMatrix<-apsOut$APSMatrix

nvar<-length(varid)
numcc<-2**nvar-1

dom<-matrix(nrow=2^nvar-1,ncol=nvar+2)
cd<-matrix(nrow=nvar-1,ncol=nvar)

varnames<-paste(rownames(PredBitMap),"-Inc",sep="")
colnames(cd)<-rownames(PredBitMap)
colnames(dom)<-c("k","R2",varnames)
rownames(dom)<-rownames(APSMatrix)
dom[,1:2]<-APSMatrix[,1:2]
 
 
for (j in 2:nvar){
  co1<-combn(varid[,1],j)
  for (k in 1:ncol(co1)){
    m<-sum(co1[,k])
    r2m<-dom[ind[m],2]
    co2<-combn(co1[,k],j-1)
      for (l in 1:ncol(co2)){
        n<-sum(co2[,l])
        r2n<-dom[ind[n],2]
        dom[ind[n],(log(m-n)/log(2))+3]<-r2m-r2n
       }
   }
}
dom<-data.frame(dom)
for (i in 1:nrow(cd)){
  ss<-subset(dom,k == i)
  cd[i,]<-colMeans(ss[,3:(nvar+2)],na.rm=TRUE)
}

cd<-insertRow(cd,1,APSMatrix[c(1:nvar),2])
rownames(cd)<-paste("CD:",c(0:(nvar-1)),sep="")
gd<-colMeans(cd)

return(list(DA=dom,CD=cd,GD=gd))
}
