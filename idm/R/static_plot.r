static_plot<-function(obj,dims=c(1,2),what=c(TRUE,TRUE),labs,pca=FALSE,contrib){
  #  require(ggplot2)
  #set to NULL to avoid compiler NOTE 
  x <- NULL
  y <- NULL
  ctr <- NULL
  
  d1 = dims[1]
  d2 = dims[2]
  out=list()
  if(pca==FALSE){
    
    if (what[2] == TRUE) {
      attdf=data.frame(x=obj$colpcoord[,d1],y=obj$colpcoord[,d2],labs=labs,ctr=obj$colctr[,d1],cor=obj$colcor[,d1],mass=obj$colmass)
      a=ggplot(data=attdf,aes(x=x,y=y))
      if (contrib == "none")
      {
        a=a+geom_text(data=attdf,aes(label=labs))+ xlab("")+ylab("") + theme(legend.position="none")
      }
      if (contrib == "cor")
      {
        a=a+geom_text(data=attdf,aes(label=labs,size=cor))+ xlab("")+ylab("")#+ theme(legend.position="none")
      }
      if (contrib == "ctr")
      {
        a=a+geom_text(data=attdf,aes(label=labs,size=ctr))+ xlab("")+ylab("")#+ theme(legend.position="none")
      }
      if (obj$ff == 0) {
        a=a+xlab(paste(round(obj$inertia_e[d1]*100,digits=2),"%",sep=""))
        a=a+ylab(paste(round(obj$inertia_e[d2]*100,digits=2),"%",sep=""))
      }
      out$attributes=a
      
    }
    
    if (what[1] == TRUE) {
      
      obsdf=data.frame(x=as.vector(obj$rowpcoord[,d1]),y=as.vector(obj$rowpcoord[,d2]),ctr=obj$rowctr[,d1],cor=obj$rowcor[,d1],mass=obj$rowmass)
      b=ggplot(data=obsdf,aes(x=x,y=y))
      b=b+geom_point(size=.5)+ xlab("")+ylab("")
      if (obj$ff == 0) { 
        b=b+xlab(paste(round(obj$inertia_e[d1]*100,digits=2),"%",sep=""))
        b=b+ylab(paste(round(obj$inertia_e[d2]*100,digits=2),"%",sep=""))
      }
      out$objects=b
      
    }
    
    if ((what[1] == TRUE) & (what[2] == TRUE)) {
      
      c=a+geom_point(data=obsdf,aes(x=x,y=y),colour="red",size=.5)+ylab("")
      if (obj$ff == 0) {
        c=c+xlab(paste(round(obj$inertia_e[d1]*100,digits=2),"%",sep=""))
        c=c+ylab(paste(round(obj$inertia_e[d2]*100,digits=2),"%",sep=""))
      }
      out$attributes=a
      out$objects=b
      out$all=c
      
    }
    
    return(out)
    
  }else{
    
    if (what[2] == TRUE) {
      attdf=data.frame(x=obj$colpcoord[,d1],y=obj$colpcoord[,d2],labs=labs)#,ctr=obj$colctr[,1],cor=obj$colcor[,1])
      cdf=circle_fun()
      a=ggplot(data=attdf,aes(x=x,y=y))
      a=a+geom_text(data=attdf,aes(label=labs))+ xlab("")+ylab("") #,size=ctr
      a=a+geom_point(data=cdf,aes(x=x,y=y),size=.05)
      a=a+geom_segment(data=attdf,aes(x=0,xend=x,y=0,yend=y),size=.5,alpha=.75)
      a=a+xlab(paste(round(obj$inertia_e[d1]*100,digits=2),"%",sep=""))
      a=a+ylab(paste(round(obj$inertia_e[d2]*100,digits=2),"%",sep=""))
      out$attributes=a 
    }
    if (what[1] == TRUE) {
      obsdf=data.frame(x=obj$rowpcoord[,d1],y=obj$rowpcoord[,d2])#,ctr=obj$rowctr[,1],cor=obj$rowcor[,1])    
      b=ggplot(data=obsdf,aes(x=x,y=y))
      b=b+geom_point(size=.5)+ xlab("")+ylab("")
      
      b=b+xlab(paste(round(obj$inertia_e[d1]*100,digits=2),"%",sep=""))
      b=b+ylab(paste(round(obj$inertia_e[d2]*100,digits=2),"%",sep=""))
      out$objects=b
    }
    
    return(out)
  }
}

circle_fun <- function(center = c(0,0),radius = 1, npoints = 1000){
  r = radius
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}