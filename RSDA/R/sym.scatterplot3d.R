sym.scatterplot3d <-
function(sym.var.x,sym.var.y,sym.var.z,labels=FALSE, ...) {
  if(((sym.var.x$var.type!='$C')||(sym.var.y$var.type!='$C')||(sym.var.z$var.type!='$C'))&&
       ((sym.var.x$var.type!='$I')||(sym.var.y$var.type!='$I')||(sym.var.y$var.type!='$I')))
    stop("Impossible to plot this type of variable")
  if((sym.var.x$var.type=='$C')&&(sym.var.y$var.type=='$C')&&(sym.var.z$var.type=='$C')) {
    if(labels==FALSE){
      p <- scatterplot3d(sym.var.x$var.data.vector,sym.var.y$var.data.vector,sym.var.z$var.data.vector, 
                         xlab=sym.var.x$var.name, ylab=sym.var.y$var.name,zlab=sym.var.z$var.name,...)
    }
    else {
      p <- scatterplot3d(sym.var.x$var.data.vector,sym.var.y$var.data.vector,sym.var.z$var.data.vector, 
                         xlab=sym.var.x$var.name, ylab=sym.var.y$var.name,zlab=sym.var.z$var.name,
                         type='n',...)
      ltext<-sym.var.x$obj.names
      text(p$xyz.convert(sym.var.x$var.data.vector,
                         sym.var.y$var.data.vector,
                         sym.var.z$var.data.vector),
           labels=ltext)
    }
  }
  if((sym.var.x$var.type=='$I')&&(sym.var.y$var.type=='$I')&&(sym.var.z$var.type=='$I')) {    
    
    xmin1<-min(sym.var.x$var.data.vector[,1])
    xmin2<-min(sym.var.x$var.data.vector[,2])
    xmin<-min(xmin1,xmin2)
    xmax1<-max(sym.var.x$var.data.vector[,1])
    xmax2<-max(sym.var.x$var.data.vector[,2])
    xmax<-max(xmax1,xmax2)
    ymin1<-min(sym.var.y$var.data.vector[,1])
    ymin2<-min(sym.var.y$var.data.vector[,2])
    ymin<-min(ymin1,ymin2)
    ymax1<-max(sym.var.y$var.data.vector[,1])    
    ymax2<-max(sym.var.y$var.data.vector[,2])    
    ymax<-max(ymax1,ymax2)
    zmin1<-min(sym.var.z$var.data.vector[,1])
    zmin2<-min(sym.var.z$var.data.vector[,2])
    zmin<-min(zmin1,zmin2)
    zmax1<-max(sym.var.z$var.data.vector[,1])    
    zmax2<-max(sym.var.z$var.data.vector[,2])    
    zmax<-max(zmax1,zmax2)
    p <- scatterplot3d(c(xmin,xmax),c(ymin,ymax),c(zmin,zmax),
                       xlab=sym.var.x$var.name, ylab=sym.var.y$var.name,zlab=sym.var.z$var.name,
                       type='n',...)
    cube <- rbind(c(1,1,1),c(2,1,1),c(2,1,2),c(1,1,2),c(1,1,1),
                  c(1,2,1),c(1,2,2),c(2,2,2),c(2,2,1),c(1,2,1),
                  c(1,2,2),c(1,1,2),c(2,1,2),c(2,2,2),c(2,2,1),c(2,1,1))
    for(i in 1:sym.var.x$N){
      
      vec.x <- sym.var.x$var.data.vector[i,cube[,1]]
      vec.y <- sym.var.y$var.data.vector[i,cube[,2]]
      vec.z <- sym.var.z$var.data.vector[i,cube[,3]]
      
      p$points3d(vec.x,vec.y,vec.z,
                 type='l',lty=1,col=i+1)
    }
    if(labels==TRUE) {
      ltext <- sym.var.x$obj.names
      textPoints <- cbind((sym.var.x$var.data.vector[,1]+sym.var.x$var.data.vector[,2])/2,
                          (sym.var.y$var.data.vector[,1]+sym.var.y$var.data.vector[,2])/2,
                          (sym.var.z$var.data.vector[,1]+sym.var.z$var.data.vector[,2])/2)
      text(p$xyz.convert(textPoints),
           labels=ltext)
    }
  }
}
