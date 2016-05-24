#' \code{indexReads} 
#' @title indexReads
#' @description for a given column of the raw count matrix this function indexes the count vector 
#'        in a vector of reads of the length of the sum of the counts. It allows after to randomly 
#'        select a subset
#' @author Edi Prifti
#' @param v : an integer raw count vector
#' @param silent : default FALSE, debugging
#' @return an character indexed vector
indexReads <- function(v, silent=FALSE){
  # keep only reads that are not 0 for optimisation purposes
  v.no0 <- v[v!=0]
  # now inflate the vector
  res <- rep(NA,sum(v.no0))
  ind.d <- 1
  ind.f <- 1
  for(i in 1:length(v.no0)){
    if(!silent){if(i%%1000==0){print(i)}}
    ind.f <- ind.d + v.no0[i]-1
    res[ind.d : ind.f] <- paste(names(v.no0)[i],1:v.no0[i],sep="_")
    ind.d <- ind.f+1
  }
  return(res)
}

#' \code{indexReadsGC} 
#' @title indexReadsGC
#' @description for a given column of the raw count matrix this function indexes the count vector 
#'        in a vector of reads of the length of the sum of the counts. It allows after to randomly 
#'        select a subset
#' @author Edi Prifti
#' @param v : an integer raw count vector
#' @param silent : default FALSE, debugging
#' @return an character indexed vector
#' @note this function is optimized when dowsizing for gene count
indexReadsGC <- function(v, silent=TRUE){
  # keep only reads that are not 0
  v.no0 <- v[v!=0]
  # now inflate the vector
  res <- rep(NA,sum(v.no0))
  ind.d <- 1
  ind.f <- 1
  for(i in 1:length(v.no0)){
    if(!silent) if(i%%1000==0){print(i)}
    ind.f <- ind.d + v.no0[i]-1
    res[ind.d : ind.f] <- names(v.no0)[i]
    ind.d <- ind.f+1
  }
  return(res)
}

#' \code{sampleRandomly} 
#' @title sampleRandomly
#' @description This function samples randomly a unique subset of the given indexed vector
#' @author Edi Prifti
#' @param v.ind : a character vector of the indexed reads
#' @param level : default 11e6, the number of reads to be selected randomly. This sould be smaller than the size of v.ind 
#' @return a character indexed vector
sampleRandomly <- function(v.ind, level=11000000){
  res <- sample(x=v.ind, size=level, replace=F)
  return(res)
}


#' \code{countSampledGenes} 
#' @title countSampledGenes
#' @description counts the number of genes that have been sampled by un_indexing
#' @author Edi Prifti
#' @param v.samp : a character vector of the sampled indexed reads (output of sampleRandlomly)
#' @return a table with counts for each gene
countSampledGenes <- function(v.samp){
  res <- gsub("_.*","",v.samp)
  res <- table(res)
  return(res)
}

#' \code{countSampledGenesGC} 
#' @title countSampledGenesGC
#' @description counts the number of genes that have been sampled by un_indexing
#' @author Edi Prifti
#' @param v.samp : a character vector of the sampled indexed reads (output of sampleRandlomly)
#' @return a table with counts for each gene
#' @note an optimized version for the gene count downsizing
countSampledGenesGC <- function(v.samp){
  res <- unique(v.samp)
  return(length(res))
}


#' \code{downsizeGC} 
#' @title downsizeGC
#' @description This function takes a matrix with raw reads counts and computes the number of genes 
#' at a given downsizing level a given number of times.
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param data : raw read count matrix with gene_ids as rownames
#' @param level : default 11e6, the downsizing level number of reads to be selected randomly. 
#' @param repetitions : default 30, the number of times the drawing is performed. Usually 30 or 10 to speed things out
#' @param silent : default is FALSE prints the status of downsizing
#' @return a matrix containing in rows a vector for each repetition and in columns the number of downsized genes for each sample
#' @note if the downsizing level is higher than the number of reads for a given sample than the result will be NA
downsizeGC <- function(data, level= 11000000, repetitions = 30, silent=FALSE){
  # for all the individuals of the set
  # data is supposed to contain raw reads and rownames to be id_fragment_external
  res <- matrix(NA,nrow=repetitions, ncol=ncol(data)); 
  rownames(res) <- 1:repetitions; colnames(res) <- colnames(data)
  for(ind in 1:ncol(data)){
    if(!silent) print(paste(ind,"Sample",colnames(data)[ind],"with",sum(data[,ind]),"reads and",sum(data[,ind]!=0),"genes"))
    if(level>sum(data[,ind])){
      if(!silent) print("This sample will not be downsized since the number of mapped reads is lower than the level")
    }else {
      v <- data[,ind]; names(v) <- rownames(data)
      v.ind <- indexReadsGC(v)
      for(step in 1:repetitions){
        v.samp <- sampleRandomly(v.ind=v.ind, level=level)
        res[step,ind] <- countSampledGenesGC(v.samp)
        if(!silent) print(paste("        step",step,"with",res[step,ind],"genes"))
      }
    }
  }
  return(res)
}


#' \code{downsizeGC.all} 
#' @title downsizeGC.all
#' @description This function takes a matrix with raw reads counts and computes the number of genes at different 
#' downsizing levels a given number of times. This is similar to the downsizeGC function but for optimization purposes 
#' it downsizes at different thresholds all together
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param data : raw read count matrix with gene_ids as rownames
#' @param levels : default seq(1E06,11E06,1E06), the downsizing levels number of reads to be selected randomly. 
#' @param repetitions : default 10, the number of times the drawing is performed. Usually 30 or 10 to speed things out
#' @param silent : default is FALSE prints the status of downsizing
#' @return a list of matrixes one per sample containing in rows a vector for each repetition and in columns the number of downsized genes for 
#'        each downsizing level
#' @note if the downsizing level is higher than the number of reads for a given sample than the result will be NA
downsizeGC.all <- function(data, levels= c(seq(1E06,11E06,1E06)), repetitions = 10, silent=FALSE){
  # for all the individuals of the set
  # data is supposed to contain raw reads and rownames to be id_fragment_external
  result <- list()
  for(ind in 1:ncol(data)){
    if(!silent) print(paste(ind,"Sample",colnames(data)[ind],"with",sum(data[,ind]),"reads and",sum(data[,ind]!=0),"genes"))
    res <- matrix(NA, nrow=repetitions, ncol=length(levels)); 
    rownames(res) <- 1:repetitions; colnames(res) <- paste("down_",levels/1E06,"M",sep="")
    for (j in 1:length(levels)){
      if(levels[j]>sum(data[,ind])){
        if(!silent) print(paste("Sample not downsized to ", levels[j]/1E06,"M since the number of mapped reads is lower",sep=""))
      }else {
        v <- data[,ind]; names(v) <- rownames(data)
        v.ind <- indexReadsGC(v)
        for(step in 1:repetitions){
          v.samp <- sampleRandomly(v.ind=v.ind, level=levels[j])
          res[step,j] <- countSampledGenesGC(v.samp)
          if(!silent) print(paste("        downsizing ",levels[j]/1E06,"M, step",step,"with",res[step,j],"genes"))
        }
      }
    }
    result[[ind]] <- res
  }
  names(result) <- colnames(data)
  return(result)
}


#' \code{downsizeMatrix} 
#' @title downsizeMatrix
#' @description takes a matrix with raw read gene counts and converts it to a downsized matrix with identical number of mapped reads 
#'      for each sample (column)
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param data : raw read count matrix with gene_ids as rownames
#' @param level : default 11E06, the downsizing levels number of reads to be selected randomly. 
#' @param repetitions : default 1, This can be also computed several times, but one is the error minimal downsizing strategy
#' @param silent : default is FALSE prints the status of downsizing
#' @return downsized read gene count matrix corresponding to the mean counts obtained with the selected number of independant 
#'      downsizing procedure
#' @note if the downsizing level is higher than the number of reads for a given sample than the result will be NA
downsizeMatrix <- function(data, level= 11000000, repetitions = 1, silent=FALSE){
  # for all the individuals of the set
  # data is supposed to contain raw reads and rownames to be id_fragment_external
  res <- matrix(NA,nrow=nrow(data),ncol=ncol(data)); rownames(res) <- rownames(data); colnames(res) <- colnames(data)
  for(ind in 1:ncol(data)){
    if(!silent) print(paste(ind,"Sample",colnames(data)[ind],"with",sum(data[,ind]),"reads and",sum(data[,ind]!=0),"genes"))
    if(sum(data[,ind]) < level){
      if(!silent) warning("This sample can't be downsized")
    }else{
      v <- data[,ind]; names(v) <- rownames(data)
      v.ind <- indexReads(v)
      # count reads for each sampled gene for each repetition and save these into a matrix
      mat.sampled <- matrix(0,nrow=nrow(data), ncol=repetitions); colnames(mat.sampled) <- 1:repetitions; rownames(mat.sampled) <- rownames(data)
      for(step in 1:repetitions){
        v.samp <- sampleRandomly(v.ind=v.ind, level=level)
        tmp <- countSampledGenes(v.samp)
        mat.sampled[match(names(tmp),rownames(mat.sampled)),step] <- tmp
        if(!silent) print(paste("        step",step,"with",sum(mat.sampled[,step]!=0),"genes"))
      }
      # compute the mean counts for each gene
      if(repetitions==1){
        res[,ind] <- mat.sampled
      }else{
        res[,ind] <- rowMeans(mat.sampled)
      }
    }
  }
  return(res)
}

#' \code{downsizedRichnessL2T} 
#' @title downsizedRichnessL2T
#' @description This procedure takes a list that is the result of the downsizeGC.all method and transforms 
#' it in a matrix of meaned downsied values. Each element of this list contains downsizing results for a given 
#' sample. This result is a matrix in lines the number of simulations and in columns the different downsizing levels
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param richness.list : the result of the downsizeGC.all method
#' @return A matrix with the samples in rows and the downsizing in columns
downsizedRichnessL2T <- function(richness.list){
  tmp <- lapply(richness.list,colMeans)
  richness.table <- c()
  for(i in 1:length(tmp)){richness.table <- cbind(richness.table ,tmp[[i]])}
  colnames(richness.table) <- names(tmp); rm(tmp)
  return(t(richness.table))
}


#' \code{computeUpsizedGC} 
#' @title computeUpsizedGC
#' @description This procedure takes a table of meaned downsized gene counts where at least one column 
#' is donwsized at a common minimal level. It uses this information to fit distributions of correlations 
#' between different downsized levels and "predict" values for the samples that have not the needed 
#' sequencing depth. The fitting of the models is based on the n-1 to be closer to reality and avoid accumulating errors.
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param richness.table : matrix with samples in rows and downsizings in the columns as produced by downsizedRichnessL2T
#' @param side : by default is 2 (downsizings in the columns)
#' @param keep.real : by default is TRUE. Substitue the predicted values by the real values when is not NA
#' @param plot : default FALSE, plots the regressions
#' @return  matrix with the same dimensions as richness.table but with complete values.
#' @seealso \code{\link{downsizedRichnessL2T}} and \code{\link{downsizeGC.all}}
computeUpsizedGC <- function(richness.table, side = 2, keep.real = TRUE, plot=FALSE){
  # regressions to be performed columnwise
  if(side!=2) richness.table <- t(richness.table)
  
  # by default we take the column having all of the points and we upsize upwords
  
  # compute pairwise regressions ans save the parameters
  richness.table <- richness.table[,order(colSums(apply(richness.table, 2, is.na)))] # order using the number of NAs
  coef <- c()
  for(i in 2:ncol(richness.table)){
    tmp <- lm(richness.table[,i] ~ richness.table[,i-1])
    if(plot) {plot(richness.table[,i] ~ richness.table[,i-1], main=colnames(richness.table)[i]); abline(tmp,col="red")}
    coef <- rbind(coef,tmp$coefficients)    
  }
  rownames(coef) <- colnames(richness.table)[-1]
  colnames(coef) <- c("intercept","slope")
  
  # compute the upsized scores
  res <- matrix(NA, nrow=nrow(richness.table), ncol=ncol(richness.table)); colnames(res) <- colnames(richness.table); rownames(res) <- rownames(richness.table)
  res[,1] <- round(richness.table[,1],0)
  for(i in 2:ncol(richness.table)){
    # There are two ways of doing this either by basing at the closest more accurate or to the first we have data
    # we chose the second option richness.table[,i]
    res[,i] <- round(res[,i-1] * coef[i-1,"slope"] + coef[i-1,"intercept"],0) # extrapolate
    
    if(keep.real){ # Subsitute the predicted values by the real ones
      res[,i][!is.na(richness.table[,i])] <- round(richness.table[,i][!is.na(richness.table[,i])],0) # keep the real data
    }
    if(plot) {plot(richness.table[,i] ~ res[,i]); abline(lm(richness.table[,i] ~ res[,i]),col="red")} # plot for checking
  }
  colnames(res) <- paste(colnames(res),"una",sep="_") # upsized NA
  return(res)
}

#' End of section and file
