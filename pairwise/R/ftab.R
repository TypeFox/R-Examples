#' @title Tabulating Answer Categories in Data
#' @export ftab
#' @description function tabulating (answer) categories in \code{X}.
#' @details \code{X} can either be a (\code{"numeric"} or \code{"character"}) \code{"matrix"} containing response vectors of persons (rows) or a \code{"data.frame"} containing \code{"numeric"}, \code{"character"} or \code{"factor"} variables (columns). 
#' @param X Data as a \code{"matrix"}, a \code{"data.frame"} or even a \code{"vector"} or \code{"factor"}. \code{"vector"} or \code{"factor"} are coerced to a \code{"data.frame"} with one column.
#' @param catgories optional a vector (\code{"numeric"} or \code{"character"}) containig the categories to tabulate. At default (\code{catgories=NULL}) the fuction looks for unique categories in \code{X}.  
#' @param na.omit logical (default: \code{na.omit=FALSE} ) wether to return frequencies for missing values, \code{NA}s.  
#' @return a \code{"matrix"} with category frequencies
#' @examples ########
#' data(bfiN)
#' ftab(bfiN)
#' data(sim200x3)
#' ftab(sim200x3)
######################################################################################

ftab <- function (X,catgories=NULL,na.omit=FALSE) {
# func. by joerg-henrik heine jhheine(at)googlemail.com  
  k<-catgories
  
  if(class(X)=="data.frame"){
    X<-sapply(X,function(x){if (class(x)=="factor"){x<-as.character(x)} else(x=x)})
    if(length(k)==0){k<-as.list(sort(unique((c(X)))))}
  }
  
  if(is.vector(X)==TRUE){X <- as.data.frame(X) }
  
  if(class(X)=="factor"){X <- as.data.frame(X) }


  # k<-as.list(sort(unique(unlist(c(X)))))}
  if(class(X)=="matrix"){
    if(length(k)==0){k<-as.list(sort(unique((c(X)))))}
  }
  # k  
  erg<-matrix(NA,ncol=dim(X)[2],nrow=length(k))
  i=1
  for (i in 1:dim(X)[2]){ 
    temp1<-mapply(rep,X[,i],length(k));dimnames(temp1)<-NULL
    sl<-as.list(data.frame(t(temp1),stringsAsFactors=F))
    # sl
    erg[,i]<-sapply(mapply(function(x,y){which(x==y)},sl,k,SIMPLIFY = FALSE),length)
  }
  rownames(erg)<-k; colnames(erg)<-colnames(X)
  if(na.omit==FALSE){miss<-apply(X,2,function(x){sum(is.na(x))})
                     erg<-rbind(miss,erg)
  }
  return(erg)
}