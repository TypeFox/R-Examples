#' @export summary.hclustgeo
#' @name summary.hclustgeo
#' @method summary hclustgeo
#' @S3method summary hclustgeo
#' @title Summary of a list of partition obtained with the method \code{hclustgeo}.
#' @description This function provide some results of several clustering obtained with \code{hclustgeo}
#' @param object an object of class 'hclustgeo' obtained with \code{hclustgeo}.
#' @param K.range the maximum number of classes we want the results.
#' @param data.desc the database used to describe the different class of the partition. If NULL the database 
#' used for the clustering will be used .
#' @param ... other arguments
#' @return {} {different results about the result of the clustering.}


summary.hclustgeo <- function(object, K.range=c(3, 4, 5), data.desc=NULL,...){
  

  #Parameters needed
  listalpha<-object[[2]]
  nb.K <- K <-length(K.range)
  nb.alpha<-length(listalpha)
  
  object<-object[[1]]
  

  res<-lapply(object, summary.hclustgeo.uniq, K.range, data.desc)
  names(res)<-paste("summary_",names(object),sep="")
  class(res)<-"summary.hclustgeo"
  
  res
}
