new.IBD<-function(snp.data, Chromosome, Position, popsam, Gen_loc=NULL, 
                  pvalue_HWE=NULL, subids=NULL,...) {
  # Set up a new object of class IBD

  # Basic checks of snp.data  
  if(!inherits(snp.data,"snp.matrix")) {
    stop("input snp.data must be of class `snp.matrix'")
  }
  if(anyDuplicated(colnames(snp.data))) {
   stop("duplicate SNP names found in input snp.data")
  }  

  # Basic checks of snp.support information
  Chromosome=gsub("[a-z]", "", Chromosome,perl=TRUE)
  snp.support<-data.frame(Chromosome=Chromosome,Position=Position)
  if (is.null(Gen_loc)){
    message("Note: Input does not include genetic map locations (Gen_loc). \nInferring genetic map from physical position (Position), \nassuming build 36 of the human genome.\n")
    snp.support$Gen_loc <- SNPgenmap(Position,Chromosome)
  } else {
    snp.support$Gen_loc <- Gen_loc
  }
  if (is.null(pvalue_HWE)) {
    message("Note: Using population sample subjects (popsam==TRUE) to fill in pvalues from tests of HWE.\n")
    sumsnp<-chopsticks::summary(snp.data[popsam,])
    snp.support$pvalue_HWE=2*pnorm(abs(sumsnp$z.HWE), lower.tail = F)
  } else {
    snp.support$pvalue_HWE<-pvalue_HWE
  }

  # Basic checks of subject.support information
  subject.support<-data.frame(popsam=popsam)
  if(is.null(subids)) {
    message("Note: Input does not include subject ids (subids). Using rownames of snp.data.\n")
    subject.support$subids<-rownames(snp.data)
  } else {
    subject.support$subids<-subids
  }

  # Ready to return an object of class IBD with our snp.data, snp.support,
  # subject.support and any information passed in via the optional argument ...
  return(IBD(snp.data,snp.support,subject.support,...))
}

IBD<-function(snp.data, snp.support,subject.support, ibd.study=NULL,
    ibd.ur=NULL,ibd.mz=NULL,ibd.po=NULL,ibd.fs=NULL,ibd.hs=NULL, 
    ibd.co=NULL, ibd.user=NULL, filterparams=NULL,simparams=NULL,call=NULL) {
  obj<-list() #Initialize an empty list, then add any ibd.* data frames we have
  obj$ibd.study=ibd.study
  obj$ibd.ur=ibd.ur
  obj$ibd.mz=ibd.mz
  obj$ibd.po=ibd.po 
  obj$ibd.fs=ibd.fs 
  obj$ibd.hs=ibd.hs
  obj$ibd.co=ibd.co
  obj$ibd.user=ibd.user
  # Now add other variables needed in the output object.
  obj<- c(obj, list(snp.data=snp.data, snp.support=snp.support,
          subject.support=subject.support,
          filterparams=filterparams, simparams=simparams, 
          call=call))
  class(obj)<-"IBD"
  return(obj)
}

is.IBD<-function(x) {
  return(inherits(x,"IBD"))
}


