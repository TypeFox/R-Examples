i_mca <- function(data1, data2=NULL,method=c("exact","live"), nchunk = 2,current_rank,ff = 0,disk=FALSE) {
  if(is.null(data2)==TRUE){
    out = mjca(data1,lambda="indicator",ret=TRUE)
    ##
    outZ = transform_z(data1,is.weight=FALSE)
    out$m = nrow(data1)
    out$rowmass = outZ$r
    out$orgn = colMeans(outZ$SZ[1:nrow(data1),])
  }else{
    if(method=="exact"){
      out=h_exact_mca(data1=data1, data2=data2,nchunk=nchunk, disk=disk)  
    }
    
    if(method=="live"){
      out=r_live_mca(data1=data1, data2=data2,nchunk=nchunk, current_rank,ff=ff,disk=disk)  
    }
  }
  class(out)="i_mca"
  return(out)
}
