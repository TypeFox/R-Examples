
.varEst<- function(H0,nFp,nPoly,pMc,pM,pMd,nCon,nCa){

pMdt<-array(NA,nPoly*nPoly)
for(i in 1:nPoly) for(j in 1:nPoly) if(i==j)pMd[i,j]=-1
count<-1
for(i in 1:nPoly) for(j in 1:nPoly) {
pMdt[count]<-pMd[i,j]
count<-count+1
}

no<-1:nFp
no.t<-no
info<-array(NA,c(choose(nPoly,nFp),nFp))
for (i in 1:choose(nPoly,nFp)){
 info[i,]<-no.t
 z<-sum(no.t==nPoly-nFp+no)
 if (z!=nFp) {
  no.t[nFp-z]<-no.t[nFp-z]+1
  if (z>0) for (j in (nFp-z+1):nFp) no.t[j]<-no.t[j-1]+1
 }
}

.C("varEstR",as.integer(info[H0,]),as.integer(nFp),as.integer(nPoly),
as.double(pMc),as.double(pM),as.double(pMdt),
as.integer(nCon),as.integer(nCa),pMcV=double(length(pMc)),PACKAGE="IFP")$pMcV
}
