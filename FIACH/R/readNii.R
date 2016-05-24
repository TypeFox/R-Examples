readNii <-function(input,fourD=TRUE){
  ################################
  ######## READ ##################
  ################################
  init<-RNiftyReg::readNifti(input)
  ##fourD TRUE,file 3d############
  
  if(length(input)>1){
    dims<-lapply(init,dim)
    ndims<-unique(sapply(dims,length))
    if(length(ndims)>1){stop("Not all images have the same number of dimensions")}
    dimMat<-do.call("rbind",dims)
    meanDim<-colMeans(dimMat)
    if(all.equal(target = meanDim,dimMat[1,])!=TRUE){stop("Not all images have the same dimensions")}
    if(length(meanDim)==4){
      if(meanDim[4]>1){stop("Will not concatenate multiple 4D images")}
    }
    if(fourD==TRUE){
      time<-length(init)
      d<-c(dimMat[1,1:3],time)
      out1<-c(init,recursive=TRUE)
      if(storage.mode(out1)=="integer"){out1<-.icombine(out1,d)}else{out1<-.dcombine(out1,d)}
      out<-RNiftyReg::updateNifti(image = out1,init[[1]])
    }else{out<-init}
  }
  
  if(length(input)==1){
        d<-dim(init)
        if(fourD==FALSE && length(d)==4){
        out<-vector(mode = "list",d[4])
        temp<-RNiftyReg::updateNifti(image = init[,,,1],template = init)
        out<-lapply(1:d[4], function(x) RNiftyReg::updateNifti(init[ , , ,x],temp))
        }else{out<-init}
        
        }
 return(out)
}



