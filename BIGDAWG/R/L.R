#' Locus Analysis Function
#'
#' This is the workhorse function for the locus level analysis.
#' @param loci.ColNames The column names of the loci being analyzed.
#' @param Locus Locus being analyzed.
#' @param genos Genotype table
#' @param grp Case/Control or Phenotype groupings.
#' @param nGrp0 Number of controls
#' @param nGrp1 Number of cases
#' @note This function is for internal BIGDAWG use only.
L <- function(loci.ColNames,Locus,genos,grp,nGrp0,nGrp1) {
  
  # pull out locus specific columns
  getCol <- seq(1,length(loci.ColNames),1)[loci.ColNames %in% Locus]
  HLA_grp <- cbind(grp,genos[,getCol])
  rownames(HLA_grp) <- NULL
  nAllele <- length(na.omit(HLA_grp[,2])) + length(na.omit(HLA_grp[,3]))
  
  ## extract alleles and counts for Grp1 and Grp0
  Alleles <- sort(unique(c(HLA_grp[,2],HLA_grp[,3])))
  
  # Build Contingency Matrix of Counts
  Allele.df <- list()      
  for(i in 1:length(Alleles)) {
    Allele.df[[i]] <- cbind(Locus,
                            Alleles[i],
                            sum(length(which(subset(HLA_grp, grp==0)[,2]==Alleles[i])),
                                length(which(subset(HLA_grp, grp==0)[,3]==Alleles[i]))),
                            sum(length(which(subset(HLA_grp, grp==1)[,2]==Alleles[i])),
                                length(which(subset(HLA_grp, grp==1)[,3]==Alleles[i]))))
  }; rm(i)
  Allele.df <- do.call(rbind,Allele.df)
  colnames(Allele.df) <- c("Locus","Allele","Group.0","Group.1")
  
  Allele.con <- matrix(as.numeric(Allele.df[,3:4]),
                       ncol=2,
                       dimnames=list(Allele.df[,'Allele'],c("Group.0", "Group.1")))
  
  if(nrow(Allele.con)>1) {
    
    ### get expected values for cells, bin small cells, and run chi square
    Result <- RunChiSq(Allele.con)
    alleles_binned <- Result$Binned
    Final_binned <- Result$Matrix
    overall.chisq <- Result$Test
    
    ## make a nice table of ORs, ci, p values
    ccdat <-TableMaker(Final_binned)
    ORout <- lapply(ccdat, cci.pval) #OR
    ORout <- do.call(rbind,ORout)
    colnames(ORout) <- c("OR","CI.lower","CI.upper","p.value","sig")
    
  } else {
    
    alleles_binned <- NA
    Final_binned <- NA
    overall.chisq <- NA
    ORout <- NA
    
  }
  
  ####################################################### Build Output List
  
  L.tmp <- list()
  
  ## Alleles.binned_out
  if(!is.null(nrow(alleles_binned))) {
    Allele.binned.tmp <- cbind(rep(Locus,nrow(alleles_binned)),
                               rownames(alleles_binned),
                               alleles_binned)
    rownames(Allele.binned.tmp) <- NULL
    colnames(Allele.binned.tmp) <- c("Locus","Allele","Group.0","Group.1")
    if(sum(grepl("\\^",Allele.binned.tmp[,'Allele']))>0) { Allele.binned.tmp[,'Allele'] <- gsub("\\^","Abs",Allele.binned.tmp[,'Allele']) }
    L.tmp[['binned']] <- Allele.binned.tmp
  } else { 
    binned.out <- cbind(Locus,'Nothing.binned',"-","-")
    colnames(binned.out) <- c("Locus","Allele","Group.0","Group.1")
    rownames(binned.out) <- NULL
    L.tmp[['binned']] <- binned.out
  }
  
  
  ## Allele.freq_out
  Allele.freq.out <- cbind(rep(Locus,nrow(Allele.con)),
                           rownames(Allele.con),
                           round(Allele.con[,'Group.0']/nGrp0,digits=5),
                           round(Allele.con[,'Group.1']/nGrp1,digits=5))
  rownames(Allele.freq.out) <- NULL
  colnames(Allele.freq.out) <- c("Locus","Allele","Group.0","Group.1")
  if(sum(grepl("\\^",Allele.freq.out[,'Allele']))>0) { Allele.freq.out[,'Allele'] <- gsub("\\^","Abs",Allele.freq.out[,'Allele']) }
  L.tmp[['freq']] <- Allele.freq.out
  
  
  ## ORtable_out
  if(!is.null(nrow(ORout))) {
    ORtable_out.tmp <- cbind(rep(Locus,nrow(ORout)),
                             rownames(ORout),
                             ORout)
    rownames(ORtable_out.tmp) <- NULL
    colnames(ORtable_out.tmp) <- c("Locus","Allele","OR","CI.lower","CI.upper","p.value","sig")
    if(sum(grepl("\\^",ORtable_out.tmp[,'Allele']))>0) { ORtable_out.tmp[,'Allele'] <- gsub("\\^","Abs",ORtable_out.tmp[,'Allele']) }
    L.tmp[['OR']] <- ORtable_out.tmp
  } else {
    ORtable_out.tmp <- cbind(Locus,Alleles,"NCalc","NCalc","NCalc","NCalc","NCalc")
    rownames(ORtable_out.tmp) <- NULL
    colnames(ORtable_out.tmp) <- c("Locus","Allele","OR","CI.lower","CI.upper","p.value","sig")
    if(sum(grepl("\\^",ORtable_out.tmp[,'Allele']))>0) { ORtable_out.tmp[,'Allele'] <- gsub("\\^","Abs",ORtable_out.tmp[,'Allele']) }
    L.tmp[['OR']] <- ORtable_out.tmp
  }
  
        
  ## overall.chisq_out
  if(!is.null(nrow(overall.chisq))) {
    overall.chisq.tmp <- cbind(Locus, overall.chisq)
    rownames(overall.chisq.tmp) <- NULL 
    colnames(overall.chisq.tmp) <- c("Locus","X.square","df","p.value","sig")
    L.tmp[['chisq']] <- overall.chisq.tmp
  } else {
    overall.chisq.tmp <- cbind(Locus,"NCalc","NCalc","NCalc","NCalc")
    rownames(overall.chisq.tmp) <- NULL 
    colnames(overall.chisq.tmp) <- c("Locus","X.square","df","p.value","sig")
    L.tmp[['chisq']] <- overall.chisq.tmp
  }
  
  
  ## Final_binned_out
  if(!is.null(nrow(Final_binned))) {
    Final_binned.tmp <- cbind(rep(Locus,nrow(Final_binned)),
                              rownames(Final_binned),
                              Final_binned)
    rownames(Final_binned.tmp) <- NULL
    colnames(Final_binned.tmp) <- c("Locus","Allele","Group.0","Group.1")
    if(sum(grepl("\\^",Final_binned.tmp[,'Allele']))>0) { Final_binned.tmp[,'Allele'] <- gsub("\\^","Abs",Final_binned.tmp[,'Allele']) }
    L.tmp[['table']] <- Final_binned.tmp
  } else {
    Final_binned.tmp <- cbind(Locus,Alleles,"NCalc","NCalc")
    rownames(Final_binned.tmp) <- NULL
    colnames(Final_binned.tmp) <- c("Locus","Allele","Group.0","Group.1")
    if(sum(grepl("\\^",Final_binned.tmp[,'Allele']))>0) { Final_binned.tmp[,'Allele'] <- gsub("\\^","Abs",Final_binned.tmp[,'Allele']) }
    L.tmp[['table']] <- Final_binned.tmp
  }
  
  
  return(L.tmp)
  

}
