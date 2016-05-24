"gplot1" <-function(M,diag=TRUE,displaylabels=TRUE,boxed.labels=FALSE,loop.cex=4,arrowhead.cex=NULL,arrowheads.fun="sqrt",edge.lwd=1,edge.col="default",rel.thresh=0.05,...){
  M[M<max(M)*rel.thresh]<-0  
  arrowheads<-t(M)[t(M)!=0]
  if(is.null(arrowhead.cex)) arrowhead.cex<-do.call(what=arrowheads.fun,list(edge.lwd))
  arrowheads<-do.call(what=arrowheads.fun,list(arrowheads))*arrowhead.cex
  if(edge.col=="default") edge.col<-gray(1-M/max(M))
  library(sna,pos=max(grep(pat="block",search()))+1)
  gplot(dat=M,diag=diag,displaylabels=displaylabels,boxed.labels=boxed.labels,loop.cex=loop.cex,arrowhead.cex=arrowheads,edge.lwd=edge.lwd,edge.col=edge.col,...)
}


"gplot2" <-
function(M,uselen=TRUE,usecurve=TRUE,edge.len=0.001,diag=TRUE,displaylabels=TRUE,boxed.labels=FALSE,loop.cex=4,arrowhead.cex=2.5,edge.lwd=1,edge.col="default",rel.thresh=0.05,...){
  M[M<max(M)*rel.thresh]<-0
  if(edge.col=="default") edge.col<-gray(1-M/max(M))
  library(sna,pos=max(grep(pat="block",search()))+1)
  gplot(dat=M,uselen=uselen,usecurve=usecurve,edge.len=edge.len,diag=diag,displaylabels=displaylabels,boxed.labels=boxed.labels,loop.cex=loop.cex,arrowhead.cex=arrowhead.cex,edge.lwd=edge.lwd,edge.col=edge.col,...)
}

