#'@title 
#'Candidate Gene Overlap
#'
#'@description 
##'\code{candidate.overlap} Extracts genes specific to samples which overlap with a candidate gene list at various pSI thresholds
#'
#'@details
#'Returns list consisting of 6 data frames, one for each pSI threshold. Each data frame contains 
#'genes specific to each sample which overlap with a candidate gene list and whose pSI values fall below each respective threshold 
#'for each cell type/sample included in the analysis.
#'NOTE:Supplementary data (human & mouse expression sets, calculated pSI datasets, etc.) can be found in \code{pSI.data} package located at the following URL:
#'\url{http://genetics.wustl.edu/jdlab/psi_package/}
#' 
#'@param pSIs  data frame output from \code{specificity.index} function with the number of columns equal to the number of samples and genes as rows.
#'@param candidate.genes candidate gene list tested for overrepresentation in cell types/samples. Comprised of official gene symbols.
#'@param write.csv logical variable indicating if csv files will be written to the current working directory (default value is FALSE)
#'
#'
#'@export
#'@import gdata
#'@author Xiaoxiao Xu, Alan B. Wells, David OBrien, Arye Nehorai, Joseph D. Dougherty
#'@examples 
#'##load sample pSI output
#'data(sample.data)
#'##load sample candidate gene lists
#'data(candidate.genes)
#'##Generates lists of overlapping genes
#'candidate.gene.overlap.AutDB <- candidate.overlap(pSIs=sample.data$pSI.output,
#'                                                  candidate.genes=candidate.genes$AutDB)
#'

candidate.overlap <- function(pSIs, candidate.genes, write.csv=FALSE){
  
  #This function utilizes the output of the pSI.list function
  list_out <- pSI.list(pSIs, write.csv=FALSE)
  
  #Range of pSI thresholds that will be used to name datasets in final output list as well as name output CSV files in the option is specified
  pSI_rang <- c("pSi_0.0001", "pSi_0.001", "pSi_0.01", "pSi_0.05")
  
  for(z in 1:length(list_out)){
    #create temporary dataset with list of genes specific to each sample for a given pSI threshold
    list_out_tmp <- list_out[[z]]
    
    #create index with match() function to identify candidate genes found also to be specific to sample at a given pSI threshold
    mtch <-function(x, y) match(y, x)
    bob <- data.frame(lapply(list_out_tmp, mtch, y=candidate.genes))
    
    #Extract lists of candidate genes found also to be specific to samples at a given pSI threshold
    f1 <-function(x, y) y[!is.na(x)]
    sxs <- data.matrix(lapply(bob, f1, y=candidate.genes))
    
    #Converts list of character arrays data frame
    for(x in 1:length(sxs)){
      if(x>1){
        pSI_holder <- cbindX(pSI_holder,data.frame(sxs[[x]]))
      }else{
        pSI_holder <- data.frame(sxs[[x]])
      }
    }
    
    #Assign column names as well as assign descriptive name to temporary dataset used in loop
    colnames(pSI_holder) <- colnames(list_out_tmp)
    assign(paste(pSI_rang[z]),pSI_holder)
    
    #Boolean used to indicate if CSV files of the output should be writted to the current working directory
    if(write.csv==TRUE){
      write.table(t(pSI_holder), paste(pSI_rang[z], "csv", sep="."), col.names=FALSE, sep=",", na="")
    }
  }
  
  #create list of 6 dataframes, one for each pSI threshold
  pSi_list <- list(get("pSi_0.0001"), get("pSi_0.001"), get("pSi_0.01"), get("pSi_0.05")) 
  names(pSi_list)<- pSI_rang
  return(pSi_list)
  
}