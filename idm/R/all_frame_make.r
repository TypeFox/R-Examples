all_frame_make<-function(obj=NULL,dims=c(1,2),nfrs=25,disk=TRUE,is.PCA=TRUE){
  outmo=list()
  for (chu in 1:obj$nchunk){
    outmo[[chu]]=frame_make(nfrs, dims=dims,disk=obj$disk,stepfra = as.numeric(chu),obj=obj,is.PCA=is.PCA)
  }
  outmo$nchunk = obj$nchunk
 # outmo$disk = obj$disk
  outmo
}



