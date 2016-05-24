#' @export plot.hclustgeo
#' @name plot.hclustgeo
#' @method plot hclustgeo
#' @S3method plot hclustgeo
#' @title Summary of a list of partition obtained with the method \code{hclustgeo}.
#' @description This function provide some results of several clustering obtained with \code{hclustgeo}
#' @param x an object of class 'hclustgeo' obtained with \code{hclustgeo}.
#' @param choice one of the following character strings: \code{'dendro'} 
#' if the user wishes to plot the dendrogram; \code{'maps'} if the user
#' wants to plot the partitions (in K.range clusters) on a map;
#' \code{'quality'} if the user wants to plot the qualities of partitions 
#' in order to choose a value for the parameter alpha.
#' @param choice.alpha a vector containig the values of parameter
#'  of alpha the user wants to plot. If \code{NULL} all the values computed in \code{hclustgeo} are used.
#' @param K.range a vector containing several values of K for which the
#'  user wishes to display plots in K clusters. If \code{NULL}, the value 
#'  of \code{K.range} is \code{c(3,4,5)}.
#' @param path.shp a character string containing the path where are the polygonial
#'  shapefiles (.shp, .dbf and .shx) to plot the maps of the typologies. If NULL, the working directory is used.
#' @param name.ind.shp a character string containing the identifier of observations in the dbf file.
#' @param save.pdf a boolean, if TRUE, graphics obtained are saved in pdf format. Default is FALSE.
#' @param path.pdf if save.pdf=TRUE, a character string containing 
#' the path where the user wishes to save generated pdf files.
#'  Default is the working directory.
#' @param ... other arguments to be passed
#' @return {} {different results and map about the result of the clustering.
#'  When the \code{plot} method is used with \code{choice = "quality"} numerical
#'   values corresponding to the quality plot are also returned.}
#' @importFrom plyr l_ply

plot.hclustgeo <-function(x, choice="maps", choice.alpha=NULL, K.range=c(3:5), 
                        path.shp=NULL, name.ind.shp=NULL,
                        save.pdf=FALSE, path.pdf=NULL, ...){
  
  
  if (!inherits(x, "hclustgeo")) 
    stop("use only with \"hclustgeo\" objects")
  
  if (!(choice %in% c("dendro", "maps", "quality"))) 
    stop("\"choice\" must be either \"dendro\", \"maps\" or \"quality\"")  
  
  res.multi<-x[[1]]
  alpha.multi<-as.double(x$multi.alpha)
  
  if (is.null(choice.alpha))
    choice.alpha <- alpha.multi
  
  #if(!all(is.element(choice.alpha, alpha.multi)))
  #  stop("\"choice.alpha\" contains values of alpha which have not been computed in \"hclustgeo\"" )
  
  
  name.res.choose <- paste0("hclustgeo.alpha=", choice.alpha) 
  res.choose<-res.multi[name.res.choose]

  
  
  ############ Plots of quality
  if (choice=="quality"){
    nb.K<-length(K.range)
    nb.alpha<-length(choice.alpha)
    
    #j is for alphas, i for clusters k, in.1 is the matrix of percentages of explained inertia 
    #(calculated on D1) according K and alpha
    #in.2 is the matrix of percentages of explained inertia calculated on D2 (geographical distance)
   
    in.1<-matrix(NA, nrow=nb.K, ncol=nb.alpha)
    rownames(in.1)<-paste("K=", K.range, sep="")
    colnames(in.1)<-paste("alpha=",choice.alpha,sep="")
    
    in.2<-in.1
    
    
    res<-lapply(res.choose, summary.hclustgeo.uniq, K.range, data.desc=NULL)
    names(res)<-paste("summary_",names(res.choose),sep="")
    
    
    for(j in 1:nb.alpha){
      for (i in 1:nb.K){
        in.1[i,j]<-res[[j]]$"heterog"[[i]]$In.var.parti
        in.2[i,j]<-res[[j]]$"heterog"[[i]]$In.geo.parti
      }
    }
    
    #Total inertia and percentages of explained inertia calculated on D1 and on D2(Dgeo) separately  
    dist.var<-res.choose[[1]]$dist.var
    dist.geo<-res.choose[[1]]$dist.geo
    wt<-res.choose[[1]]$wt
    
    Itot.var<-heterog.Ck(dist.var, dist.geo, partition=rep(1,length(wt)), classe=1, wt, alpha=1)$In.var.Ck
    Itot.geo<-heterog.Ck(dist.var, dist.geo, partition=rep(1,length(wt)), classe=1, wt, alpha=1)$In.geo.Ck
    
    
    Itot.1<-Itot.var
    Itot.2<-Itot.geo
    
    pourc.inertie.1<-1-(in.1/Itot.1)
    pourc.inertie.2<-1-(in.2/Itot.2)
    
    mycol<-1:nb.K
    myx<-0:(nb.alpha-1)
    
    
    #Plot of percentages of explained inerta calculated on D1
    #We also return numerical values on a data frame
    if(save.pdf==TRUE){
      setwd(path.pdf)
      pdf(file="Percentages of explained inertia on variables.pdf",width=8.26,height=11.7)
    }
    
    quali.part.var<-matrix(rep(NA,nb.K*nb.alpha),nrow=nb.K,ncol=nb.alpha)
    rownames(quali.part.var)<-paste0("K=",K.range)
    colnames(quali.part.var)<-paste0("alpha=",choice.alpha)
    
    
    par(mfrow=c(1,1))
    plot(0,0,xaxt="n",ylim=c(0,1),xlim=c(0,(nb.alpha-1)),type="n",xlab="alpha",
         ylab="Percentage of explained inertia",main="Percentages of explained inertia on variables")
    axis(1, at=myx, labels=choice.alpha)
    for(i in 1:nb.K){
      quali.part.var[i, ]<-pourc.inertie.1[i,]
      points(x=myx,y=pourc.inertie.1[i,],col=mycol[i],pch=16,type="b")
    }
    legend("topleft",legend=rownames(pourc.inertie.1), fill=mycol)
    
    
    if(save.pdf==TRUE){dev.off()}
    
    
    #Plot of percentages of explained inerta calculated on D2
    #We also return numerical values on a data frame
    if(save.pdf==TRUE){
      setwd(path.pdf)
      pdf(file="Percentages of explained inertia on geographical distances.pdf",width=8.26,height=11.7)
    }
    
    quali.part.geo<-matrix(rep(NA,nb.K*nb.alpha),nrow=nb.K,ncol=nb.alpha)
    rownames(quali.part.geo)<-paste0("K=",K.range)
    colnames(quali.part.geo)<-paste0("alpha=",choice.alpha)
    
    
    par(mfrow=c(1,1))
    plot(0,0,xaxt="n",ylim=c(0,1),xlim=c(0,(nb.alpha-1)),type="n",xlab="alpha",
         ylab="Percentage of explained inertia",main="Percentages of explained inertia on geographical distances")
    axis(1, at=myx, labels=choice.alpha)
    for(i in 1:nb.K){
      quali.part.geo[i, ]<-pourc.inertie.2[i,]
      points(x=myx,y=pourc.inertie.2[i,],col=mycol[i],pch=16,type="b")
    }
    legend("topright",legend=rownames(pourc.inertie.2), fill=mycol)
    if(save.pdf==TRUE){dev.off()}
    
    quali.part<-list(quali.part.var=quali.part.var, quali.part.geo=quali.part.geo)
    return(quali.part)
    
  }
  ############ End of Plots of quality
  
  
  ############ Plots of dendrograms
  if (choice=="dendro"){
    
    lapply(res.choose, plot.hclustgeo.uniq, choice="dendro", K.range=K.range, path.shp= path.shp,
           name.ind.shp=name.ind.shp, save.pdf=save.pdf, path.pdf=path.pdf)
    
  }
  ############ End of Plots of dendrograms
  
  ############ Plots of maps
  if (choice=="maps"){
    
    lapply(res.choose, plot.hclustgeo.uniq, choice="maps", K.range=K.range, path.shp= path.shp,
           name.ind.shp=name.ind.shp, save.pdf=save.pdf, path.pdf=path.pdf)
    
  }
  ############ End of Plots of maps


}
