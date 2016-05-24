#' \code{normFreqRPKM} 
#' @title normFreqRPKM
#' @description converts a raw count matrix onto a frequency matrix using the RPKM normalization method.
#'        This method consists of two consecutive steps, first dividing the raw counts by the length of 
#'        the gene sequence and the second shrinking the signal per column to a sum of 1
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param dat : raw counts data matrix with gene_ids as rownames
#' @param cat : the current working catalogue where the reads are mapped and counted, (i.e. hs_3.3_metahit, hs_3.9_metahit)
#'        This can also be a vector of genelength values that correspond to the number of rows in the dat matrix and are 
#'        ordered respectively
#' @return a normalized frequency matrix
normFreqRPKM <- function(dat, cat = NULL){
  # load the annotation datafile to extract the fragment size
  if(is.null(cat)){
    stop("cat should be provided as a vector containing the gene length of the catalog with the same names as the profile matrix.")
  }else{
    if(length(cat) != nrow(dat)){
      stop("cat should have the same length as the number of rows in the dat profile matrix.")
    }else{
      genesize <- cat[match(rownames(dat),names(cat))]
    }
  }
  
  if(is.matrix(dat)) print("The dataset is a matrix")
  else{ # transform the data.frame onto a matrix
    dat <- as.matrix(dat)
  }
  # divide by the genesize
  res <- dat
  for(i in 1:ncol(dat)) {
    res[,i] <- dat[,i]/genesize
  }  
  # divide by the total number of reads
  for(i in 1:ncol(res)) {
    res[,i] <- res[,i]/sum(res[,i])
  }
  return(res)
}


#' \code{normFreqTC} 
#' @title normFreqTC
#' @description converts a raw count matrix onto a frequency matrix using the TC (total count) normalization method.
#'        This method consists of scaling the signal by the total counts per each sample
#' @author Edi Prifti
#' @param dat : raw counts data matrix with gene_ids as rownames
#' @return a normalized frequency matrix
normFreqTC <- function(dat){
  if(is.matrix(dat)) print("The dataset is a matrix")
  else{ # transform the data.frame onto a matrix
    dat <- as.matrix(dat)
  }
  # divide by the genesize
  res <- dat
  # divide by the total number of reads
  for(i in 1:ncol(res)) {
    res[,i] <- res[,i]/sum(res[,i])
  }
  return(res)
}

#' ADD other normalization algorithms
#' End of section and file