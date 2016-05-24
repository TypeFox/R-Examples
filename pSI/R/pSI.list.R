#'@title 
#'Convert pSI output to gene list
#'
#'@description 
#'\code{pSI.list} returns list consisting of 6 data frames, one for each pSI threshold.
#'
#'@details
#'Each data frame contains genes whose pSI values fall below each respective threshold for each cell type/sample included in the analysis.
#'NOTE:Supplementary data (human & mouse expression sets, calculated pSI datasets, etc.) can be found in \code{pSI.data} package located at the following URL:
#'\url{http://genetics.wustl.edu/jdlab/psi_package/}
#'
#'@param pSIs  data frame output from \code{specificity.index} function with the number of columns equal to the number of samples and genes as rows.
#'@param write.csv logical variable indicating if csv files will be written to the current working directory (default value is FALSE)
#'
#'@export
#'@import gdata
#'@author Xiaoxiao Xu, Alan B. Wells, David OBrien, Arye Nehorai, Joseph D. Dougherty
#'
#'@examples
#'##load sample pSI output
#'data(sample.data)
#'##List the genes specific to each cell type/sample type across all pSI thresholds
#'pSI.out.list <- pSI.list(pSIs=sample.data$pSI.output, write.csv=FALSE)
#'

pSI.list <- function(pSIs,write.csv=TRUE){
  
  #Reassign gene names as uppercase for continuity in downstream analysis
  rownames(pSIs) <- toupper(rownames(pSIs))
  
  #Assigning variables used in loop, ie number of genes, number of samples etc.
  gene_l <- length(pSIs[,2])
  n_samples <- length(pSIs[1,])
  
  
  for(b in c(0.05, 0.01, 0.001, 0.0001)){
    
    for(i in 1:n_samples){
      #Remove genes with pSI values of NA, ie ignore genes whose SI values did not fall within the top 10% of all genes
      sig_tmp <- pSIs[!is.na(pSIs[,i]),]
      #Keep genes with pSI values below threshold, ie ignore genes whose pSI values were above set pSI threshold value
      sig_tmp <- sig_tmp[sig_tmp[,i] < b,]
      sig_tmp <- data.frame(rownames(sig_tmp), stringsAsFactors=FALSE)
      
      #Binds together the lists of genes specific to each sample type for a given threshold
      if(!i==1){
        gene_bind <- cbindX(gene_bind, sig_tmp)
        
      }     
      else{
        gene_bind <- sig_tmp
      }  
    }
    #Assigns columns names to temporary loop data frame as well assign temporary loop data frame a descriptive name
    #colnames(gene_bind) <- paste(colnames(pSIs), b, sep="_")
    if(b<0.001){
      colnames(gene_bind) <- paste(colnames(pSIs), "0.0001", sep="_")
    }else{
      colnames(gene_bind) <- paste(colnames(pSIs), b, sep="_")
    }
    if(b<0.001){
      assign(paste("pSi", "0.0001", sep="_"),gene_bind)
    }else{
      assign(paste("pSi", b, sep="_"),gene_bind)
    }
    #assign(paste("pSi", b, sep="_"),gene_bind)
    
    #Boolean used to indicate if CSV files of the output should be writted to the current working directory
    if(write.csv==TRUE){
      write.table(t(gene_bind), paste("pSi", paste(b,".csv",sep=""), sep="_"), col.names=FALSE, sep=",", na="")
    }
    
    if(b<0.05){
      out <- c(out,gene_bind) 
    }else{
      out <-gene_bind
    }  
    
  }
  
  #create list of 6 dataframes, one for each pSI threshold
  pSi_list <- list(get("pSi_0.0001"), get("pSi_0.001"), get("pSi_0.01"), get("pSi_0.05"))
  pSI_rang <- c("pSi_0.0001", "pSi_0.001", "pSi_0.01", "pSi_0.05")
  names(pSi_list)<- pSI_rang
  
  return(pSi_list)

}


