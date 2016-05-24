#' Haplotype Analysis Function
#'
#' This is the workhorse function for the haplotype analysis.
#' @param genos.sub The genotype columns of the loci(locus) set being analyzed. 
#' @param grp Case/Control or Phenotype groupings.
#' @note This function is for internal BIGDAWG use only.
H <- function(genos.sub,grp) {

    loci.sub <- unique(gsub(".1","",colnames(genos.sub),fixed=T))
    nloci.sub <- as.numeric(length(loci.sub))
    Haplotype <- paste(loci.sub,collapse="~")  
    
    cat("Estimating Haplotypes ...",Haplotype,"\n")
    
    ### estimate haplotypes
    Tab.out <- haplo.stats::haplo.em(geno=genos.sub,locus.label=loci.sub)
    
    ## extract haplotype freqs for cases and controls
    Subjects <- as.list(seq_len(nrow(genos.sub)))
    Tab.Haps <- lapply(Subjects,FUN=getHap,HaploEM=Tab.out)
    Tab.Haps <- cbind(grp,do.call(rbind,Tab.Haps))
    colnames(Tab.Haps)[2:3] <-c("Haplotype.1","Haplotype.2")
    
    ## Build Contingency Matrix of Counts
    haps <- sort(unique(c(Tab.Haps[,'Haplotype.1'],Tab.Haps[,'Haplotype.2'])))
    haps_counts <- mat.or.vec(nr=length(haps),nc=2)
    rownames(haps_counts) <- haps
    colnames(haps_counts) <- c('Group.0','Group.1')
    
      #Group 0
      Tab.Haps.grp0  <- Tab.Haps[which(Tab.Haps[,'grp']==0),c('Haplotype.1','Haplotype.2')]
      haps_grp0 <- table(Tab.Haps.grp0)
      haps_counts[match(names(haps_grp0),haps),'Group.0'] <- haps_grp0
      
      #Group 1
      Tab.Haps.grp1  <- Tab.Haps[which(Tab.Haps[,'grp']==1),c('Haplotype.1','Haplotype.2')]
      haps_grp1 <- table(Tab.Haps.grp1)
      haps_counts[match(names(haps_grp1),haps),'Group.1'] <- haps_grp1
    
    ### get expected values for cells, bin small cells, and run chi square
    Result <- RunChiSq(haps_counts)
    haps_binned <- Result$Binned
    Final_binned <- Result$Matrix
    overall.chisq <- Result$Test
   
    ## Convert counts to frequencies
    haps_freq <- haps_counts
    haps_freq[,'Group.0'] <- haps_freq[,'Group.0']/(nrow(Tab.Haps.grp0)*2)
    haps_freq[,'Group.1'] <- haps_freq[,'Group.1']/(nrow(Tab.Haps.grp1)*2)
     
    ## make a nice table of ORs, ci, p values
    ccdat <-TableMaker(Final_binned)
    ORout <- lapply(ccdat, cci.pval) #OR
    ORout <- do.call(rbind,ORout)
    
    #haps_binned - Binned Haplotypes
    haps_binned_fix <- cbind(rownames(haps_binned),haps_binned)
    colnames(haps_binned_fix) <- c(Haplotype,colnames(haps_binned))
    rownames(haps_binned_fix) <- NULL
    if(sum(grepl("\\^",haps_binned_fix[,Haplotype]))>0) { haps_binned_fix[,Haplotype] <- gsub("\\^","Abs",haps_binned_fix[,Haplotype]) }
    
    #final_binned - Contingency Table for ChiSq
    Final_binned_fix <- cbind(rownames(Final_binned),Final_binned)
    colnames(Final_binned_fix) <- c(Haplotype,colnames(Final_binned))
    rownames(Final_binned_fix) <- NULL
    if(sum(grepl("^",Final_binned_fix[,Haplotype]))>0) { Final_binned_fix[,Haplotype] <- gsub("\\^","Abs",Final_binned_fix[,Haplotype]) }
    
    #haps_freq - Frequencies
    haps_freq_fix <- cbind(rownames(haps_freq),haps_freq)
    colnames(haps_freq_fix) <- c(Haplotype,colnames(haps_freq))
    rownames(haps_freq_fix) <- NULL
    if(sum(grepl("\\^",haps_freq_fix[,Haplotype]))>0) { haps_freq_fix[,Haplotype] <- gsub("\\^","Abs",haps_freq_fix[,Haplotype]) }
    
    #ORout - ODDs Ratios
    ORout_fix <- cbind(rownames(ORout),ORout)
    colnames(ORout_fix) <- c(Haplotype,colnames(ORout))
    rownames(ORout_fix) <- NULL
    if(sum(grepl("\\^",ORout_fix[,Haplotype]))>0) { ORout_fix[,Haplotype] <- gsub("\\^","Abs",ORout_fix[,Haplotype]) }
    
    #Haplotype - Replace Abs symbols
    if(sum(grepl("\\^",Tab.Haps[,'Haplotype.1'])) + sum(grepl("\\^",Tab.Haps[,'Haplotype.2'])) >0) { 
      Tab.Haps[,'Haplotype.1'] <- gsub("\\^","Abs",Tab.Haps[,'Haplotype.1'])
      Tab.Haps[,'Haplotype.2'] <- gsub("\\^","Abs",Tab.Haps[,'Haplotype.2'])
      }
    colnames(Tab.Haps)[2:3] <- c(paste(Haplotype,".Hap1",sep=""),paste(Haplotype,".Hap2",sep=""))
    
    H.tmp <- list()
    H.tmp[['Haplotypes']] <- Tab.Haps[,2:3] # table of subject haplotypes
    H.tmp[['freq']] <- haps_freq_fix # rounded frequencies  
    H.tmp[['binned']] <- haps_binned_fix # binned haplotypes
    H.tmp[['OR']] <- ORout_fix # odd ratio table
    H.tmp[['chisq']] <- overall.chisq # chi sq test statistic
    H.tmp[['table']] <- Final_binned_fix # final table for chi sq

    return(H.tmp)
  
}
