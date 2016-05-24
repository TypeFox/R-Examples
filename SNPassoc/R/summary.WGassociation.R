`summary.WGassociation` <-
function(object,pSig=0.000001,...)
{

 if (!inherits(object, "WGassociation"))
      stop("object must be of class 'WGassociation'")

 genes<-attr(object,"gen.info")

 if (is.null(genes)){
   aux<-attr(object,"pvalues")
   info0<-nrow(aux)
   info1<-round((table(factor(aux[,1],levels=c("Genot error","Monomorphic")))/info0)*100,1)
   nSig<-sum(aux[,2]<=pSig,na.rm=TRUE)
   info2<-c(nSig,round((nSig/info0)*100,1))
   info<-c(info0,info1,info2)
   ans<-rbind(info)
   rownames(ans)<-""
   colnames(ans)<-c("SNPs (n)","Genot error (%)","Monomorphic (%)","Significant* (n)","(%)")
} 
else {
  SNPs<-attr(object,"label.SNPs")
  pval<-attr(object,"pvalues")
  nSNPs<-table(genes[,2])
  chr.l <- names(nSNPs)

  if(length(chr.l)==22)
   o<-orderChromosome(chr.l)
  else 
   o<-c(1:length(chr.l))
  chr <- chr.l[o]
  nSNPs.o<-nSNPs[o]

  info<-matrix(NA,nrow=length(chr),ncol=5)

  for (i in 1:length(chr))
    {
     info0<-nSNPs.o[i]
     temp<-genes[genes[,2]==chr[i],]
     aux<-pval[dimnames(pval)[[1]]%in%temp[,1],]
     info1<-round((table(factor(aux[,1],levels=c("Genot error","Monomorphic")))/nrow(aux))*100,1)
     nSig<-sum(aux[,2]<=pSig,na.rm=TRUE)
     info2<-c(nSig,round((nSig/nrow(aux))*100,1))
     info[i,]<-c(info0,info1,info2)
    }

  ans<-data.frame(info)
  names(ans)<-c("SNPs (n)","Genot error (%)","Monomorphic (%)",
      "Significant* (n)","(%)")
  dimnames(ans)[[1]]<-chr

}
print(ans)
cat("\n *Number of statistically significant associations at level", pSig)
cat("\n")
invisible(ans)
}

