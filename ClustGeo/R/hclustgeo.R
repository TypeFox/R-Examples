#' @export
#' @name hclustgeo
#' @title \code{hclustgeo.uniq} for several parameters alpha.
#' @description Perform \code{hclustgeo.uniq} for several parameter alpha.
#' @param data a data frame with \code{n} rows and \code{p} colums where \code{n} observations are
#'  described by \code{p} numerical variables.
#' @param D.geo a matrix of size \code{n} by \code{n} containing geographical distances between observations
#' @param alpha a vector containing all the values of the parameter alpha for which the user wishes to perform \code{hclustgeo}.
#' @return {}{an object of class hclustgeo containing all the results of \code{hclustgeo.uniq} for the differents values of \code{alpha}}
#' @examples
#' ###load data
#' #library("ClustGeo")
#' #data(comm303)
#' #base <- comm303$data.303
#' #Dgeo <- comm303$Dgeo.303
#' 
#' ###perform hclustgeo for alpha=1
#' #res.alpha1 <- hclustgeo(data=base, D.geo=Dgeo, alpha=1)
#' 
#' ###plot of the dendrogram
#' #plot(res.alpha1, choice="dendro") #we retain the partition in K=5 clusters
#' 
#' ###plot of the map
#' #path.303 <- file.path(path.package("ClustGeo"), "shapes/comm303")
#' #ID.ind <- "INSEE_COM"
#' #plot(res.alpha1, choice="maps", K.range=5, path.shp=path.303, name.ind.shp=ID.ind)
#' 
#' ###perform hclustgeo for several values of alpha
#' #multi.alpha <- seq(0, 1, 0.2)
#' #res.alpha <- hclustgeo(data=base, D.geo=Dgeo, alpha=multi.alpha)
#' 
#' ###plot of the qualities
#' #plot.qual <- plot(res.alpha, choice="quality", K.range=5,
#' #                    path.shp=path.303, name.ind.shp=ID.ind)




#Fonction pour boucler sur tous les alphas
hclustgeo <- function(data, D.geo, alpha){
  multi.alpha <- as.vector(alpha) #alpha is a vector we call it multi.alpha in the function
  multi.alpha<-sort(multi.alpha)
  nb.alpha<-length(multi.alpha)
  list.res<-list()
  for(i in 1: nb.alpha){
    cat("#### hclustgeo_alpha=",multi.alpha[i],"\n",sep="")
    list.res[[i]]<-hclustgeo.uniq(data, D.geo, alpha=multi.alpha[i])
  }
  names(list.res)<-paste("hclustgeo.alpha=",multi.alpha,sep="")
  
  list.res2<-list(results.hclustgeo=list.res, multi.alpha=multi.alpha)
  class(list.res2)<-"hclustgeo"
  return(list.res2)
}
