computePenalty<-function(readSupport, readWeights, pUnknown){
totalReads<-sum(readWeights[,"weight"])
pij<-matrix(data=0, ncol=2, nrow=totalReads)
LLunknown<-totalReads*log(pUnknown)
pij[,1]<-pUnknown
pij[1:readSupport,2]<-1
abund<-c((totalReads - readSupport)/totalReads, readSupport/totalReads)
lpenalty<-LLunknown-colSums(log(pij%*%abund))
return(lpenalty)
}                      

computePenalty.nucl<-function(readSupport, readWeights, pUnknown, median.genome.length){
totalReads<-sum(readWeights[,"weight"])
pij<-matrix(data=0, ncol=2, nrow=totalReads)
LLunknown<-totalReads*log(pUnknown)
pij[,1]<-pUnknown
pij[1:readSupport,2]<-1/median.genome.length
abund<-c((totalReads - readSupport)/totalReads, readSupport/totalReads)
lpenalty<-LLunknown-colSums(log(pij%*%abund))
return(lpenalty)
}                      
