plot.i_mca<-function(x,dims=c(1,2),what=c(TRUE,TRUE),contrib="none",dataname=NULL,
                     labels=NULL,animation=TRUE,frames=10,zoom=TRUE, ...){
  
  # require(animation)
  #  source("static_plot.R")
  #  source("all_frame_make.R")
  #  source("ani_plot.R")
  obj <- x
  #default dataname
  if (is.null(dataname)){
    dataname = deparse(substitute(x))
  }
  
  #default labels
  if (is.null(labels)){
    labels = obj$levelnames
  }

  if(animation==TRUE){
    outmo=all_frame_make(obj=obj,dims=dims,nfrs=frames,is.PCA=FALSE)  
    
    #attributes
    if (what[2] == TRUE) {
      movieNameA=paste("iMCA",dataname,"atts","movie.gif",sep="_")
      ani_plot(moname=movieNameA,outmo,nfrs=frames,labs=labels,att=TRUE,pca=FALSE,contrib=contrib,zoom=zoom)
    }
    
    #objects
    if (what[1] == TRUE) {
      movieNameO=paste("iMCA",dataname,"obs","movie.gif",sep="_")
      ani_plot(moname=movieNameO,outmo,nfrs=frames,labs=labels,att=FALSE,pca=FALSE,contrib=contrib,zoom=zoom)
    }
  }else{
    stpl=static_plot(obj, dims=dims, what=what,labs=labels,pca=FALSE,contrib=contrib)
    return(stpl)
  }
}

