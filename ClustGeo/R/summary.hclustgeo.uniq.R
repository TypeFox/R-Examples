#' @export summary.hclustgeo.uniq
#' @name summary.hclustgeo.uniq
#' @method summary hclustgeo.uniq
#' @S3method summary hclustgeo.uniq
#' @title Summary of a partition obtained with the method \code{hclustgeo.uniq}.
#' @description This function provide some results of a clustering obtained with \code{hclustgeo.uniq}
#' @param object an object of class 'hclustgeo.uniq' obtained with \code{hclustgeo.uniq}.
#' @param K.range the maximum number of classes we want the results.
#' @param data.desc the database used to describe the different class of the partition. If NULL the database 
#' used for the clustering will be used .
#' @param ... other arguments
#' @return {} {different results about the result of the clustering.}
#' @keywords internal
#' @importFrom  FactoMineR catdes



summary.hclustgeo.uniq <- function(object,  K.range=c(3:5), data.desc=NULL,...){
  
  res<-object
  
  if(is.null(data.desc)){
    base<-res$data
  }else{base<-data.desc}
  
  
  dist.var<-res$dist.var
  dist.geo<-res$dist.geo
  wt<-res$wt
  
  
  
  alpha<-res$alpha
  nom<-rownames(base)
  p<-ncol(base)
  nb.K<-length(K.range)
  
  
  #List of partitions in K.range clusters
  list.cut<-list()
  for (i in 1: nb.K){
    list.cut[[i]]<-cutree(res,K.range[i])
  }
  names(list.cut)<-paste("K=",K.range,sep="")
  
  
  #List of quality of partitions in K.range clusters
  list.heterog<-list()
  for (i in 1: nb.K){
    list.heterog[[i]]<-heterog.parti(dist.var, dist.geo, partition=list.cut[[i]], wt, alpha=alpha)
  }
  names(list.heterog)<-paste("K=",K.range,sep="")
  
  
  
  #Description of classes of the different partitions
  base.typo<-cbind(nom,data.frame(list.cut))
  list.desc<-list()
  for(i in 1:nb.K){
    base.desc<-data.frame(base,as.factor(base.typo[,i+1]))
    list.desc[[i]]<-catdes(base.desc,p+1)$quanti
  }
  names(list.desc)<-paste("K=",K.range,sep="")
  
  
  
  resultat<-list(alpha=alpha, K.range=K.range, cut=base.typo, heterog=list.heterog, desc=list.desc, data.desc=base)
  class(resultat) <- "summary.hclustgeo.uniq"
#   
#   toprint <- matrix("",7,2)
#   colnames(toprint) <-c("name","description")
#   toprint[1,] <- c("$alpha", "parameter 'alpha' used in the ClustGeo method")
#   toprint[2,] <- c("$K.range", "parameter 'K.range' used in the method") 
#   toprint[3,] <- c("$cut", "dataframe whose colums correspond to the different typologies in 'K.range' clusters")
#   toprint[4,] <- c("$heterog", "list of heterogeneity criterions for each partitions")
#   toprint[5,] <- c("$desc", "List of descriptions of clusters for each partition generated.")
#   toprint[6,] <- c("","Partitions are described by variables of 'data.desc'")  
#   toprint[7,] <- c("$data.desc", "The data base used to describe cluster of partitions")    
#  
#   print(toprint)
#   
  
  return(resultat)
  
}
