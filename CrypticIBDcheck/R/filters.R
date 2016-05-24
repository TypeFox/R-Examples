## filter.control is a function to set parameters that control quality control
## filters: snpcallrate, MAF, samplecallrate and HWEp 

filter.control<-function(filter=TRUE,snpcallrate=0.9,MAF=0.01,
                         samplecallrate=0.9, HWEp=0.001) {
  list(filter=filter,snpcallrate=snpcallrate,MAF=MAF,
       samplecallrate=samplecallrate, HWEp=HWEp)
}

################################################################
snpfilter <- function(snpmatlist1,filter){

SNP.support=snpmatlist1$snp.support
SNP.support=remove.cdlIBS(SNP.support)  #remove cdlIBS columns, if present
snpobjt=snpmatlist1$snp.data

sumsnp<-chopsticks::summary(snpobjt)
nbs=ncol(snpobjt)

#snps call rate test
 callr=(sumsnp$Call.rate>=filter$snpcallrate) 
 snpobjt1<-snpobjt[,(callr)]
 SNP.support1<-SNP.support[callr,]

#MAF test
 sumsnp1<-chopsticks::summary(snpobjt1)
 maft=(sumsnp1$MAF>=filter$MAF) 
 snpobjt1<-snpobjt1[,(maft)]
 SNP.support1<-SNP.support1[maft,]

#samples call rate test
 callrvt=(row.summary(snpobjt1)$Call.rate>=filter$samplecallrate) 
 snpobjt1<-snpobjt1[(callrvt),]
 snpmatlist1$subject.support<-snpmatlist1$subject.support[(callrvt),]

#Remove SNPs not found by SNPgenmap (NA genetic location)
  ind<-is.na(SNP.support1$Gen_loc)
  SNP.support1<-SNP.support1[!ind,]
  snpobjt1<-snpobjt1[,(!ind)]
   
#Do test of HWE and keep SNPs that pass
  vecpval=as.numeric(as.vector(SNP.support1$pvalue_HWE))
  hwvec<-( vecpval>filter$HWEp & !is.na(SNP.support1$pvalue_HWE) )
   
  snpmatlist1$snp.data<-snpobjt1[,(hwvec)]
  snpmatlist1$snp.support<-SNP.support1[hwvec,]
  
 
    return(snpmatlist1)
    
}

