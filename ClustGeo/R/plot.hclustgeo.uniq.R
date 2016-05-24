#' @export plot.hclustgeo.uniq
#' @name plot.hclustgeo.uniq
#' @method plot hclustgeo.uniq
#' @S3method plot hclustgeo.uniq
#' @title Differents plots available with the method \code{hclustgeo.uniq}.
#' @description This function provide some graphics and maps obtained with the method \code{hclustgeo.uniq}
#' @param x an object of class 'hclustgeo.uniq' obtained with \code{hclustgeo.uniq}.
#' @param choice a character string: "dendro" if you want to plot the dendrogram; "maps" if you want to plot
#' the partitions (in K.range clusters) on a map.
#' @param K.range a vector containing the different partitions (each K.range correspond to a number of cluster)
#' you want to represent on a map.
#' @param path.shp the path where are the shp and dbf files to plot the map of the typologie.
#' @param name.ind.shp the identifier of individuals in the dbf file.
#' @param save.pdf a boolean indicating if you want to save the maps obtained in pdf files.
#' @param path.pdf if save.pdf=TRUE, the path of the folder where you want to save pdf files.
#' @param ... other arguments to be passed
#' @return {} {different results and map about the result of the clustering.}
#' @keywords internal
#' @importFrom rCarto mapChoroTypo


plot.hclustgeo.uniq <- function(x, choice="maps", K.range=c(3:5), 
                             path.shp=NULL, name.ind.shp=NULL,
                             save.pdf=FALSE, path.pdf=NULL,...){
  
  #if (!inherits(x, "hclustgeo.uniq")) 
  # stop("use only with \"hclustgeo.uniq\" objects")
  
  if (!(choice %in% c("dendro", "maps"))) 
    stop("\"choice\" must be either \"dendro\" or \"maps\"")  
  
  if (choice=="maps" & (is.null(path.shp)|is.null(name.ind.shp)))
    stop("If you want to plot maps, you have to provide path to shp files in \"path.shp\" and identifiers of individuals in \"name.ind.shp\".")
  
  if (save.pdf==TRUE & is.null(path.pdf))
    path.pdf<-getwd()
  
  res<-x
  

  base<-res$data

  
  alpha<-res$alpha
  nom<-rownames(base)
  p<-ncol(base)
  nb.K<-length(K.range)
  

  if (choice=="dendro"){
    
    #Plot of dendrogram of res
    res.hclust<-list(merge=res$merge, height=res$height, order=res$order)
    class(res.hclust)<-"hclust"
    
    title.plot<-paste("Dendrogram_hclustgeo_alpha=",alpha,sep="")
    title.pdf<-paste("Dendrogram_hclustgeo_alpha=",alpha,".pdf",sep="")
    if(save.pdf==TRUE){
      setwd(path.pdf)
      pdf(file=title.pdf,width=8.26,height=11.7)
    }
    par(mfrow=c(1,1))
    dendro<-plot(res.hclust,main=title.plot)
    dendro.plot<-recordPlot()
    if(save.pdf==TRUE){dev.off()}
    
    list.map <-NULL
    
  }
  


  if (choice=="maps"){
    
    #List of partitions in K.range clusters
    list.cut<-list()
    for (i in 1: nb.K){
      list.cut[[i]]<-cutree(res,K.range[i])
    }
    names(list.cut)<-paste("K=",K.range,sep="")
    
    
    #Plot of maps in K.range clusters
    list.map<-list()
    base.plot<-cbind(nom,data.frame(list.cut))
    
    for (i in 1:nb.K){
      K.choose<-K.range[i]
      
      title.plot<-paste("Map_hclustgeo_alpha=",alpha,"\n",K.choose,"_clusters",sep="")
      title.pdf<-paste("Map_hclustgeo_alpha=",alpha,"_",K.choose,"clusters.pdf",sep="")
      
      if(save.pdf==TRUE){
        adresse.new<-paste(path.pdf,"/",K.choose,"_classes",sep="")
        dir.create(adresse.new, showWarnings=FALSE)
        setwd(adresse.new)
        pdf(file=title.pdf,width=8.26,height=11.7)
      }
      
      #Here we want to supresse some warnings of mapChorotypo
      #specially warning for plots in 2 clusters
      h <- function(w) if( any( grepl( "minimal value for n is 3", w) ) ) invokeRestart( "muffleWarning" ) 
      
      withCallingHandlers(
      mapChoroTypo(shpFile=path.shp, shpId=name.ind.shp,
                   df=base.plot,dfId="nom",palCol="Paired",var=colnames(base.plot)[i+1],width=15,
                   height=20,title=title.plot,...), warning=h)
      
      
      list.map[[i]]<-recordPlot()
      if(save.pdf==TRUE){dev.off()}
    }
    
    names(list.map)<-paste("Map_",K.range,"_clusters",sep="")
    dendro.plot <- NULL
    
  }

  resultat<-list(alpha=alpha, maps=list.map, dendro=dendro.plot, K.range=K.range,
                 path.shp=path.shp, path.pdf=path.pdf, name.ind.shp=name.ind.shp)



  class(resultat) <- "plot.hclustgeo.uniq"
  resultat
  
}
