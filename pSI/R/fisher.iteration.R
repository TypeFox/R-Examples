#'@title 
#'Fisher's Exact Test Across All Cell Types & pSI Thresholds
#'
#'@description 
#'\code{fisher.iteration} will test a candidate gene list for overrepresenation in the various cell type/pSI threshold combinations produced by
#'the specificty.index function. 
#'NOTE:Supplementary data (human & mouse expression sets, calculated pSI datasets, etc.) can be found in \code{pSI.data} package located at the following URL:
#'\url{http://genetics.wustl.edu/jdlab/psi_package/}
#'
#'@details
#'This function is used to answer the question of what is the probability that a certain number of genes specific to a certain cell type/sample occured by chance
#'(as usual with low probabilities corresponding to high statistical significance).
#'This is accomplished with a binary variable for each gene in the population with two mutual exclusive values: 1) The gene is specific to the 
#'cell type/sample in question or 2) The gene is not specific to the cell type/sample in question
#'
#'@param pSIs  data frame output from \code{specificity.index} function with the number of columns equal to the number of samples and genes as rows.
#'@param candidate.genes candidate gene list tested for overrepresentation in cell types/samples. Comprised of official gene symbols.
#'@param background character string used to indicate what background gene list should be used in Fisher's exact test for overrepresentation. 
#'The default value is \code{"data.set"} which indicates that the gene list of the input pSI data set will be used to represent the background gene list. 
#'This would be used in the case when the input pSI data set is comprised of genes derived from the same species as the genes found in the candidate gene list.
#'\code{background} can take on two other values, the first of which is \code{"human.mouse"}. 
#'\code{"human.mouse"} indicates that the background gene list will be comprised of intersection of two lists: 1) all genes in the input pSI dataset (all are human genes), 2) all genes with clear human-mouse homologs.
#'This option would be used in the case when the input data set is comprised of human genes (i.e. genes from a human microarray) and the candidate gene list being tested is comprised of mouse genes. 
#'The last value \code{background} can take on is  \code{"mouse.human"}.
#'\code{"mouse.human"} indicates that the background gene list will be comprised of intersection of two lists: 1) all genes in the input pSI dataset (all are mouse genes), 2) all genes with clear mouse-human homologs.
#'This option would be used in the case when the input data set is comprised of mouse genes (i.e. genes from a mouse microarray) and the candidate gene list being tested is comprised of human genes. 
#'@param p.adjust logical. default output is bonferroni corrected p-value but if \code{p.adjust} is \code{FALSE}, nominal p-values will be output.
#'
#'@export
#'@author Xiaoxiao Xu, Alan B. Wells, David OBrien, Arye Nehorai, Joseph D. Dougherty
#'@examples
#'##load sample pSI output
#'data(sample.data)

#'##load sample candidate gene lists
#'data(candidate.genes)
#'##run Fisher's exact test for overrperesentation on pSI.out for the AutDB 
#'##candidate gene list across all cell types/sample types & pSI thresholds
#'fisher.out.AutDB <- fisher.iteration(pSIs=sample.data$pSI.output, 
#'                                          candidate.genes=candidate.genes$AutDB)
#'

fisher.iteration<-function(pSIs, candidate.genes, background="data.set", p.adjust=TRUE){  
  
  if(background=="data.set"){
    
    total <- length(pSIs[,1])
    
  }else if(background=="human.mouse"){
    
    temp <- tempfile()
    download.file("http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/kgXref.txt.gz", temp)
    kgxref.mouse <- read.table(gzfile(temp), sep="\t", row.names=NULL, header=FALSE, stringsAsFactors=FALSE, quote="")
    unlink(temp)
    kgxref.mouse <- kgxref.mouse[,5]
    total <- sum(!is.na(match(toupper(rownames(pSIs)), toupper(kgxref.mouse))))
    
  }else if(background=="mouse.human"){
    
    temp <- tempfile()
    download.file("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz", temp)
    kgxref.human <- read.table(gzfile(temp), sep="\t", row.names=NULL, header=FALSE, stringsAsFactors=FALSE, quote="")
    unlink(temp)
    kgxref.human <- kgxref.human[,5]
    #total <- 13421
    total <- sum(!is.na(match(toupper(rownames(pSIs)), toupper(kgxref.human))))
    
  }else{
    stop("Background gene list for Fisher's exact test incorrectly assigned")
  }
  
  #Assigning variables used in loop
  pSIs<-cbind(c(1:length(pSIs[,1])),pSIs)
  index <- cbind(as.numeric(c(1:length(pSIs[,1]))),toupper(rownames(pSIs)))
  
  #create 3 column dataset which is comprised of gene index, boolean if gene is found in candidate list, and gene symbol in the third column
  candidate.genes <- cbind(index[,1], match(index[,2], toupper(unique(candidate.genes))), index[,2])
  
  #Create temporary variable, cnv_temp and assign value of zero to genes with no match and value of 1 to genes with a match
  cnv_temp <- candidate.genes
  cnv_temp[which(is.na(candidate.genes[,2])),2] <- 0
  cnv_temp[which(!is.na(candidate.genes[,2])),2] <- 1
  
  #remove index of genes from pSI dataset and assign rownames to index(NULL)
  pSIs <- pSIs[,-1]
  row.names(pSIs) <- NULL
  
  #reassign cnv_temp variable as candidate.genes while dropping the third column and drop temporary variable
  candidate.genes <- cnv_temp[,-3]
  rm(cnv_temp)
  
  #range of pSI thresholds to loop through
  rang<- c(.05,.01,.001,.0001)
  
  #matrix created hold output from loop
  outs_CTp<-data.frame(matrix(NA, nrow=length(colnames(pSIs)), ncol=length(rang)))
  rownames(outs_CTp) <- colnames(pSIs)
  colnames(outs_CTp) <- rang
  
  #Loops through the various samples/pSI threshold combinations and test for overrepresentation using a Fisher's Exact test
  j<-0
  for (Th in rang){
    j<-j+1
    for(i in 1:length(colnames(pSIs))){
      outs_CTp[i,j]<-as.numeric(fisher(pSIs, candidate.genes, Th, i, total))
    }
  }
  if(p.adjust==TRUE){
  
  #Create new dataframe "holder" which has twice as many columns as number of samples/cell types.
  colnames(outs_CTp) <- paste(colnames(outs_CTp),"adjusted",sep=" - ")
  
  #Calculte BH adjusted p-values and add to "holder" dataframe
  for(k in 1:ncol(outs_CTp)){
    outs_CTp[,k] <- p.adjust(outs_CTp[,k], method="BH", n=nrow(outs_CTp))
  }

  }else{
    colnames(outs_CTp) <- paste(colnames(outs_CTp),"nominal",sep=" - ")
  }
  
  return(outs_CTp)
}