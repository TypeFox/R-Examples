plot.i_pca<-function(x,dims=c(1,2),what=c(TRUE,TRUE),dataname=NULL,
                     labels=NULL, animation=TRUE,frames=10,zoom=TRUE, ...){
  
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
    outmo=all_frame_make(obj=obj,dims=dims,nfrs=frames,is.PCA=TRUE)
    
    #attributes
    if (what[2] == TRUE) {
      movieNameA=paste("iPCA",dataname,"atts","movie.gif",sep="_")
      ani_plot(moname=movieNameA,outmo,nfrs=frames,labs=labels,att=TRUE,pca=TRUE,zoom=zoom)
    }
    
    #objects
    if (what[1] == TRUE) {
      movieNameO=paste("iPCA",dataname,"obs","movie.gif",sep="_")
      ani_plot(moname=movieNameO,outmo,nfrs=frames,labs=labels,att=FALSE,pca=TRUE,zoom=zoom)
    }
  }else{
    stpl=static_plot(obj,dims=dims,what=what,labs=labels,pca=TRUE)
    return(stpl)
  }
  
  
}