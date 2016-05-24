###########################################################################
#####                           SegCorr                               #####
###########################################################################
### Summary: Performs block diagonal segmentation. There is an option to correct the gene expression signal 
###          for CNV. Then, segmentation and the exact test are performed.

### Input:   CHR:         an id vector stating to which chromosome each gene belongs to. The vector should have the same
###                       length as the number of rows of the EXP matrix.
###          EXP:         (Raw) Gene expression matrix. The number of rows should account for the genes and the number of columns 
###                       for the patients
###          genes:       the name of genes.The vector should have the same
###                       length as the number of rows of the EXP matrix.
###          Position.EXP: a matrix that should have the same number of rows as the EXP matrix and have 2 columns. The 
###                        first column is the start postion and the second is the end position.
###          CNV:          is a logical variable. if CNV=T then CNV correction is performed.
###          SNP.CHR:      (variable provided only when CNV=T) an id vector stating to which chromosome each SNP
###                        probe belongs to. The vector should have the same
###                        length as the number of rows of the SNP matrix.
###          SNP:          (variable provided only when CNV=T) SNP matrix. The number of rows should account for the SNP and the number of columns 
###                       for the patients. The ordering of the patients must be the same as in the EXP matrix.
###         Position.SNP:  matrix that should have the same length as the as the SNP probes in the SNP  matrix.
###         group:         variance classification of the SNP probes. if variable missing, all patients are assumed to be
###                         allocated in the same group.

### Output: Region.List:   a matrix containing the highly correlated regions in descending order (the most significant on top).
###                        CHR: the chromosome of the region. Start/End: are the order of the gene within each chromosome.
###                        Rho: is the correlation.  length: number of genes in the region.
###                        first/last gene: the name of the first/last gene in the region
###                        pvalue: the pvalue as obtained from the test. genes: the names of the genes belonging to the region
###                        pvalue.adj: the pvalue corrected for multiple testing.
###         
###       Chromosome.Inf:  a matrix containning the estimated background correlation (rho0.hat) per chromsome, the number of segments
###                        and the loglikehood.

#options(warn=1)
SegCorr  = function(CHR,EXP,genes,S,CNV,SNPSMOOTH,Position.EXP,SNP.CHR,SNP,Position.SNP,group){
  
  
  
  if(missing(EXP)){   
    stop("Expression matrix missing")     
  }else{
    p = dim(EXP)[1]
    n = dim(EXP)[2]
  }
  
  if(missing(CHR)){
    stop('Chromosome vector missing') 
  }else{
    if(length(CHR)!=p){
      stop("Dimension mismatch: Expression matrix and chromosomes")
    }
  }
  
  if(missing(genes)){
    warning('Gene names missing: using artificial names')
    genes = 1:p
  }else{
    if(length(genes)!=p){
      stop("Dimension mismatch: Expression matrix and gene names")
    }
  }
  
  if(missing(CNV)){
    warning('No CNV correction performed: Default value used')
    CNV = F
  }else if(!missing(CNV) & !is.logical(CNV)){
    stop('CNV: Input must be True/False')
  }
  
  if(missing(S)){
    S=0.7
    warning("Using default value for S")
  }
  
  
  if(CNV == F){
    
    chromosome = unique(CHR)
    Region.List = vector('list',length(chromosome))
    Chromosome.Inf = matrix(NA,length(chromosome),3)
    for(i in 1:length(chromosome)){
      loc = which(CHR==chromosome[i])
      x = segmentation(unique(CHR[loc]),EXP[loc,],genes[loc])
      Region.List[[i]] = x$Results
      Chromosome.Inf[i,] = c(x$rho0,x$K,x$likelihood) 
    }
    
    colnames(Chromosome.Inf) = c('rho0.hat',"no of segments","likelihood")
    rownames(Chromosome.Inf) = chromosome
    
    Region.List = do.call("rbind",Region.List)
    Region.List$p.values.adj = p.adjust(Region.List$p.value,method ='BH')
    Region.List = Region.List[order(Region.List$p.values.adj),]
    rownames(Region.List) = NULL
    list(Region.List = Region.List, Chromosome.Inf = Chromosome.Inf)
    
  }else if (CNV == T){
    
    if(missing(Position.EXP)){
      stop('Position.EXP missing')
    }else{
      if(dim(Position.EXP)[1]!=p){
        stop("Dimension mismatch: Expression matrix and Position matrix")
      }else if(dim(Position.EXP)[2]!=2){
        stop("Position.EXP: must be a matrix with 2 columns")
      }
    }
    
    
    if(missing(SNP)){   
      stop("SNP matrix missing")     
    }else{
      p.SNP = dim(SNP)[1]
      n.SNP = dim(SNP)[2]
      
      if(n.SNP !=n ){
        stop('Patient mistmatch between EXP and SNP')
      }
    }
    
    
    if(missing(SNP.CHR)){
      stop('SNP.CHR vector missing') 
    }else{
      if(length(SNP.CHR)!=p.SNP){
        stop("Dimension mismatch: SNP matrix and SNP.CHR")
      }else if(length(intersect(unique(SNP.CHR),unique(CHR)) )!= length(unique(CHR))){
        stop('SNP matrix does not contain information for all Chromosomes needed')
      }
    }
    
    
    if(missing(Position.SNP)){
      stop('Position.SNP missing')
    }else{
      if(length(Position.SNP)!=p.SNP){
        stop("Dimension mismatch: SNP matrix and Position matrix")
      }
    }
    
    if(missing(SNPSMOOTH)){
      warning('No SNP mean segmentation performed: Default value used')
      SNPSMOOTH = F
    }else if(!missing(SNPSMOOTH) & !is.logical(SNPSMOOTH)){
      stop('SNPSMOOTH: Input must be True/False')
    }
    
    
    
    
    
    
    
    chromosome = unique(CHR)
    Region.List = EXP.corrected = vector('list',length(chromosome))
    Chromosome.Inf = matrix(NA,length(chromosome),3)
    for(i in 1:length(chromosome)){
      cat(paste('Now analysing Chromosome =',chromosome[i],sep=" "),sep="\n")
      loc = which(CHR==chromosome[i])
      loc.SNP = which(SNP.CHR==chromosome[i])
      
      if(SNPSMOOTH==T){
        
        if(missing(group)){
          group = rep(1,n)
        }else{
          if(length(group)!=n){
            stop("Dimension mismatch: group and number of patients")
          }
        }
        
        mu.SNP = segmented_signal(SNP[loc.SNP,] ,group)
      }else{
        cat(paste('No SNP mean segmentation performed'),sep="\n")
        mu.SNP = SNP[loc.SNP,]
      }
      
      
    
      correction = CNV_correction(Position.EXP[loc,1],Position.EXP[loc,2], Position.SNP[loc.SNP], mu.SNP, EXP[loc,])
      
      EXP.corrected[[i]] = correction
      x = segmentation(unique(CHR[loc]),correction,genes[loc])
      Region.List[[i]] = x$Results
      Chromosome.Inf[i,] = c(x$rho0,x$K,x$likelihood) 
    }
    
    colnames(Chromosome.Inf) = c('rho0.hat',"no of segments","likelihood")
    rownames(Chromosome.Inf) = chromosome
    EXP.corrected = do.call("rbind",EXP.corrected)
    Region.List = do.call("rbind",Region.List)
    
    Region.List$p.values.adj = p.adjust(Region.List$p.value,method ='BH')
    Region.List = Region.List[order(Region.List$p.values.adj),]
    rownames(Region.List) = NULL
    list(Region.List = Region.List, Chromosome.Inf = Chromosome.Inf, EXP.corrected = EXP.corrected )
    
    
    
    
    
  }
  
  
  
  
  
}
