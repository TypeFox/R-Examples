#' Creates a sample
#'
#' R Implementation of the SPSS \code{SAMPLE} argument. Takes a sample from a xpssFrame object, data frame or matrix.
#'
#' xpssSample takes a sample of the specified size from the elements of x either with or without replacement. The subset get specified by pct or n. 
#' 
#' \code{pct} specifies a percentage value, for the amount of data which should be kept, allowed value range is from 0.01 to 1. \cr
#' \code{n} indicates the amount of values to keep. \cr
#' \code{from} determines the basis for n. \code{from} has to be higher then \code{n}. \cr
#'
#' @param x a (non-empty) data.frame or input data of class \code{"xpssFrame"}. 
#' @param pct atomic numeric, determines the percentage to keep.
#' @param n atomic numeric, specifies the number of cases to keep.
#' @param from atomic numeric, indicates the Basis for n.
#' @return Returns a subset of the actual dataset. 
#' @author Andreas Wygrabek
#' @seealso \code{\link{sample}}
#' @examples
#' data(fromXPSS)
#' 
#' xpssSample(fromXPSS, pct = 0.5)
#' @export
xpssSample <- function(x, pct = NULL, n = NULL, from = NULL){
    
    ####################################################################
    ####################################################################
    
    functiontype <- "DM"
    x <- applyMetaCheck(x)
    
    ####################################################################
    ####################################################################
    ####################################################################
    
    
    
    attr_backup <- attributesBackup(x)
    
    if(!is.null(pct) && !is.null(n))
    {
      stop("it is only possible to use pct or n, but not both arguments at the same time.")
    }
    if(!is.null(from)){
      if(length(x[[1]]) < from)
      {
        stop("from is longer than the dataset, chose a basis value which is smaller or has a equal length of the dataset")
      } 
      if(is.null(n)) 
      {
        stop("if from is set, n has to be set, too")
      }
    }
    if(!is.null(n)){
      if(length(x[[1]]) <= n)
      {
        stop("n is longer than the dataset, chose a basis value which has a smaller length of the dataset")
      }
      if(is.null(from)) 
      {
        stop("if n is set, from has to be set, too")
      }
    }
    if(!is.null(pct)) {
      if((pct < 0) || (pct>1))
      {
        stop("pct value range is limited from 0 to 1")
      }  
    }
    if(!is.null(n) && !is.null(from)) {
      if(n > from)
      {
        stop(" n < from is not true")
      }
    }
     
      if(!is.null(pct)){
        if(is.data.frame(x) | is.matrix(x)){
          
          sampVec <- 1:nrow(x)
          
          sampVec <- sample(sampVec,round(length(sampVec)*pct, 0))
          
          out <- x[sampVec,]
          
        } else {
          
          sampVec <- 1:length(x)
          
          sampVec <- sample(sampVec,round(length(sampVec)*pct, 0))
          
          out <- x[sampVec]
          
        }
      } else if(is.null(from)) {
        if(is.data.frame(x) | is.matrix(x)){
          
          sampVec <- 1:nrow(x)
          
          sampVec <- sample(sampVec,n, 0)
          
          out <- x[sampVec,]
          
        } else {
          
          sampVec <- 1:length(x)
          
          sampVec <- sample(sampVec,n, 0)
          
          out <- x[sampVec]
          
        }   
      } else if(!is.null(n) & !is.null(from) & if(is.data.frame(x) | is.matrix(x)){nrow(x)>from} else {length(x)>from}){
        if(is.data.frame(x) | is.matrix(x)){
          
          sampVec <- 1:nrow(x)
          
          sampVec <- sample(sampVec[c(1:from)],round(length(sampVec[c(1:from)])*(n/from),0), 0)
          
          out <- x[sampVec,]
          
          } else {
          
          sampVec <- 1:length(x)
          
          sampVec <- sample(sampVec[c(1:from)],round(length(sampVec[c(1:from)])*(n/from),0), 0)
          
          out <- x[sampVec]
          
        } 
        
      }
      else {
        if(is.data.frame(x) | is.matrix(x)){
          
          sampVec <- 1:nrow(x)
          
          sampVec <- sample(sampVec,round(length(sampVec)*(n/from),0), 0)
          
          out <- x[sampVec,]
          
        } else {
          
          sampVec <- 1:length(x)
          
          sampVec <- sample(sampVec,round(length(sampVec)*(n/from),0), 0)
          
          out <- x[sampVec]
        }   
      }  
    
    for(i in 1:length(out)) {
      attributes(out[[i]]) <- attr_backup$local[[i]]  
  
      na_pos <- which(rownames(out[i]) %in% attributes(out[[i]])$MIS)
      value_pos <- which(attributes(out[[i]])$MIS %in% rownames(out[i]))
      attr_temp <- attributes(out[[i]])$MIS
      
      if(length(na_pos>0)){ 
        attributes(out[[i]])$MIS <- NULL
        attributes(out[[i]])$MIS <- data.frame("POS" = na_pos, "VAL"=attr_temp[,2][value_pos])
      }
      
    } 
    
    
    
    if(length(attributes(out)$FILTER)>0 && attributes(out)$FILTER != FALSE) {
      pos <- which(!(eval(parse(text = paste("out$",(attributes(x)$FILTER),sep="")))))
      attributes(out)$FILTERED_DATA <- out[pos,]
      pos <- which(eval(parse(text = paste("out$",(attributes(x)$FILTER),sep=""))))
      out <-  out[pos,]      
    } else { 
      x <- applyAttributeDemerge(x)
    }

    return(out)
  }
    
   
