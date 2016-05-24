
#' \code{filterListGenebags} 
#' @title filterListGenebags
#' @description This function filters a list of mgs with gene_id >= a given mgs gene number.
#' @author Emmanuelle Le Chatelier
#' @param list.genebags : a list of genebags that needs to be filtered
#' @param size.min : the minimal size threshold above which mgs are selected, default 0
#' @param size.max : the maximal size threshold above which mgs are selected, default 15000
#' @return a list of selected genebags with their original elements (usually geneids)
#' @note This is the former filterListMGS function
filterListGenebags <- function(list.genebags, size.min = 0, size.max = 15000){ 
  if(!is.list(list.genebags)){
    stop("Error! this is not a list")
  }
  res <- list.genebags[as.numeric(summary(list.genebags)[,1]) >= size.min & as.numeric(summary(list.genebags)[,1]) <= size.max ]
  return(res)
}


#' \code{filterMat} 
#' @title filterMat
#' @description This function filters a matrix (mat) based on the rate of positive value (filt) 
#' in each individual under a given degree of presence (filt), all sparse values are put to 0.
#' @author Edi Prifti & Emmanuelle Le Chatelier
#' @param mat : matrix of counts 
#' @param filt : filtering threshold in percentage 
#' @return a cleaned matrix
filterMat <- function(mat, filt=0){
  if(!is.matrix(mat)){
    stop("Error! only a matrix can be provided ")
  }
  tmp <- mat
  tmp[,colSums(mat!=0)<round(nrow(mat)*filt/100)] <- 0
  return(tmp)
}

#' End of section and file