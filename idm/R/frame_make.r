frame_make<-function(nfrs,dims=c(1,2),stepfra=1,disk=TRUE,obj=NULL,is.PCA=TRUE){
  a = dims[1]
  b = dims[2]
  if(disk==TRUE){
    fStart=paste("./PCstory/PCstart",deparse(stepfra),".txt",sep="")
    fEnd=paste("./PCstory/PCEnd",deparse(stepfra+1),".txt",sep="")
    ufStart=paste("./PCstory/PCstartUnit",deparse(stepfra),".txt",sep="")
    ufEnd=paste("./PCstory/PCEndUnit",deparse(stepfra+1),".txt",sep="")

    if(is.PCA!=TRUE){
      snam.colctr=paste("./PCstory/PCctr",deparse(stepfra),".txt",sep="")
      enam.colctr=paste("./PCstory/PCctr",deparse(stepfra+1),".txt",sep="")
      snam.colcor=paste("./PCstory/PCcor",deparse(stepfra),".txt",sep="")
      enam.colcor=paste("./PCstory/PCcor",deparse(stepfra+1),".txt",sep="")
    }
    
    snam.nrowctr=paste("./PCstory/PCctrUnit",deparse(stepfra),".txt",sep="")
    enam.nrowctr=paste("./PCstory/PCctrUnit",deparse(stepfra+1),".txt",sep="")
    snam.rowcor=paste("./PCstory/PCcorUnit",deparse(stepfra),".txt",sep="")
    enam.rowcor=paste("./PCstory/PCcorUnit",deparse(stepfra+1),".txt",sep="")
    
    dat1=read.table(fStart)
    dat2=read.table(fEnd)
    #  udat1=read.table(ufStart,row.names=NULL)[,c(a+1,b+1)]
    #  udat2=read.table(ufEnd,row.names=NULL)[,c(a+1,b+1)]
    udat1=read.table(ufStart,row.names=NULL)[,-1]
    udat2=read.table(ufEnd,row.names=NULL)[,-1]
    newrow=nrow(udat2)-nrow(udat1)
    nudat=matrix(0,newrow,dim(udat1)[2])
    udat1=rbind(udat1,nudat)
    if(is.PCA!=TRUE){
      if (length(snam.colctr) != 0) {
        scolctr=read.table(snam.colctr)
        ecolctr=read.table(enam.colctr)
      }  
      if (length(snam.colcor) != 0) {
        scolcor=read.table(snam.colcor)
        ecolcor=read.table(enam.colcor)
      }
    }
    srowctr=read.table(snam.nrowctr,row.names=NULL)[,-1]
    erowctr=read.table(enam.nrowctr,row.names=NULL)[,-1]
    
    srowcor=read.table(snam.rowcor,row.names=NULL)[,-1]
    erowcor=read.table(enam.rowcor,row.names=NULL)[,-1]
    
  }else{
    dat1=obj$allcolcoord[[stepfra]]
    dat2=obj$allcolcoord[[stepfra+1]]
    udat1=obj$allrowcoord[[stepfra]]
    udat2=obj$allrowcoord[[stepfra+1]]
    newrow=nrow(udat2)-nrow(udat1)
    nudat=matrix(0,newrow,ncol(dat1))
    udat1=rbind(udat1,nudat)
    
    if(is.PCA!=TRUE){
      # absolute contributions
      scolctr=obj$allcolctr[[stepfra]] 
      ecolctr=obj$allcolctr[[stepfra+1]]
      
      # relative contributions
      scolcor=obj$allcolcor[[stepfra]]
      ecolcor=obj$allcolcor[[stepfra+1]]
    }
  }
  
  if(is.PCA!=TRUE){
    if (length(scolctr) != 0) {
      colctr=cbind(((sqrt(scolctr[,a])+sqrt(scolctr[,b]))^2)/2,((sqrt(ecolctr[,a])+sqrt(ecolctr[,b]))^2)/2)
    }
    if (length(scolcor) != 0) {
      colcor=cbind(((sqrt(scolcor[,a])+sqrt(scolcor[,b]))^2)/2,((sqrt(ecolcor[,a])+sqrt(ecolcor[,b]))^2)/2)
    }
  }
  #print(colctr)
  #dat=cbind(dat1[,1],dat1[,2],dat2[,1],dat2[,2])
  dat=cbind(dat1[,a],dat2[,a],dat1[,b],dat2[,b])
  udat=cbind(udat1[,a],udat2[,a],udat1[,b],udat2[,b])
  xran=range(dat[,1:2])
  yran=range(dat[,3:4])
  uxran=range(udat[,1:2])
  uyran=range(udat[,3:4])
  outmo=list()
  outmo$dat=morph(dat,frames = nfrs)
  outmo$udat=morph(udat,frames = nfrs)
  outmo$xr=xran
  outmo$yr=yran
  outmo$uxr=uxran
  outmo$uyr=uyran  
  
  if(is.PCA!=TRUE){
    if (length(scolctr) != 0) {
      outmo$bfcolctr=apply(colctr,1,function(x) (seq(from=x[1],to=x[2],length=nfrs)))
    }
    if (length(scolcor) != 0) {
      outmo$bfcolcor=apply(colcor,1,function(x) (seq(from=x[1],to=x[2],length=nfrs)))
    }
  }
  
  outmo
}
