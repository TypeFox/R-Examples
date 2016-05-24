getProbsRegions.i<-function(i,blocks,probs,annotation,nclass)
 {
   probes<-dimnames(probs)[[1]]
   cond<-annotation$kb>=blocks[i,2] & annotation$kb<=blocks[i,3] & annotation$Chrom==blocks[i,1]
   selec<-annotation[cond,1]
   ans<-probs[probes%in%selec,]
   p<-ans[1,-c(1:4)]
   out<-matrix(as.numeric(p),ncol=nclass,byrow=TRUE)
   out
 }

