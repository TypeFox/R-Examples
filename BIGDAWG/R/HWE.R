#' Hardy Weinbergy Equilibrium Function
#'
#' This is the main wrapper function for each analysis.
#' @param Tab data frame of genotype files post processing.
#' @param All.ColNames character vector of Tab object column names.
#' @note This function is for internal BIGDAWG use only.
HWEChiSq <- function(Tab,All.ColNames) {
  
  genos.sub <- Tab[which(Tab[,2]==0),3:ncol(Tab)]
  loci <- as.list(unique(All.ColNames[3:length(All.ColNames)]))
  nloci <- length(loci)
  
  #Format genotypes
  df.1 <- data.frame(genos.sub[,seq(1,nloci*2,2)])
  df.2 <- data.frame(genos.sub[,seq(2,nloci*2,2)])
  colnames(df.2) <- colnames(df.1)
  df <- rbind(df.1,df.2)
  colnames(df) <- do.call(rbind,loci)
  rm(df.1,df.2)
  
  #Allele info
  Alleles <- lapply(loci,FUN=function(x) sort(unique(df[,x]))); names(Alleles) <- loci # unique allele names
  nAlleles <- lapply(loci,FUN=function(x) length(na.omit(unique(df[,x])))); names(nAlleles) <- loci # no. unique alleles
  nAlleles.tot <- lapply(loci,FUN=function(x) length(df[,x])); names(nAlleles.tot) <- loci # total no. alleles
  
  #Possible Genotypes
  Allele.Combn <- lapply(nAlleles,makeComb); names(Allele.Combn) <- loci
  
  #Get Allele Counts and Frequencies
  Allele.cnts <- lapply(loci,FUN=function(x) table(df[,x])); names(Allele.cnts) <- loci
  Allele.Freq <- lapply(loci,FUN=function(x) Allele.cnts[[x]]/nAlleles.tot[[x]]); names(Allele.Freq) <- loci
  
  #Get Observed and Expected Frequencies Matrix for genotypes
  Freq.Final <- lapply(loci,FUN=getCS.Mat,genos.sub=genos.sub,Allele.Freq=Allele.Freq,Allele.Combn=Allele.Combn)
  names(Freq.Final) <- loci
  
  #Calculate Chi Square Statistic for each Locus
  Freq.chisq <- lapply(loci,FUN=getCS.stat,Freq.Final=Freq.Final)
  names(Freq.chisq) <- loci
  
  #Recompute number of alleles and genotypes at each locus from binned contigency matrices (Freq.Final)
  nAlleles.bin <- lapply(Freq.Final,FUN=getAllele.Count)
    names(nAlleles.bin) <- loci
  nGenotypes.bin <- lapply(Freq.Final,nrow)
    names(nGenotypes.bin) <- loci
    nGenotypes.bin[which(as.numeric(lapply(nGenotypes.bin,FUN=is.null))==1)] <- NA
  
  #Get degrees of freedom for each locus
  #Alleles = a ; possible genotypes =g ; df = g - (a - 1)
  Allele.dof <- lapply(loci,FUN=function(x) nGenotypes.bin[[x]] - (nAlleles.bin[[x]] - 1) ); names(Allele.dof) <- loci
  
  #Get P.values from Chi Square distribution
  Freq.pvals <- lapply(loci,FUN=function(x) 1-pchisq(as.numeric(Freq.chisq[[x]]), as.numeric(Allele.dof[[x]])))
  names(Freq.pvals) <- loci
  
  #Format Output
  Test.out <- cbind(loci,
                    do.call(rbind,Freq.chisq),
                    do.call(rbind,Allele.dof),
                    do.call(rbind,Freq.pvals),
                    rep('NS',nloci))
  colnames(Test.out) <- c("Locus","X.square","df","p.value","sig")
  rownames(Test.out) <- NULL
  
  Test.out[which(as.numeric(Test.out[,'p.value'])<0.05),"sig"] <- "*"
  Test.out <- matrix(unlist(Test.out),ncol=ncol(Test.out),dimnames=dimnames(Test.out))
  rownames(Test.out) <- NULL
  
  Test.out[,'X.square'] <- sapply(as.numeric(Test.out[,'X.square']),FUN=round,digits=4)
  Test.out[,'p.value'] <- sapply(as.numeric(Test.out[,'p.value']),FUN=function(x) format.pval(x))

  #Flag for invalid degrees of freedom
  flagLoci <- which(as.numeric(Test.out[,'df'])<1)
    
  #Flag for invalid chi square matrices
  Freq.Flag <- lapply(loci,FUN=function(x) ifelse(nrow(Freq.Final[[x]])>2,0,1))
  flagLoci <- unique(c(flagLoci,which(Freq.Flag==1)))
  
  #Flag for invalid chi square matrices
  flagLoci <- unique(c(flagLoci,which(is.na(Test.out[,'X.square']))))
  
  if(length(flagLoci)>0){
    Test.out[Test.out[,'Locus'] %in% unlist(loci[flagLoci]),2:ncol(Test.out)] <- "NCalc"
  }  
  
  return(Test.out)
  
}
