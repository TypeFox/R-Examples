delim.xtab<-function(x,pct=c("row","column","cell"),coltot=TRUE,rowtot=TRUE,
 ndec=1,delim="\t",interdigitate=TRUE,label=deparse(substitute(x))) {

 xmat<-x$counts
 if(coltot) {
  xmat<-rbind(x$counts,x$col.margin)
  rownames(xmat)<-c(rownames(x$counts),"Total")
 }
 if(rowtot) {
  xmat<-cbind(xmat,c(x$row.margin,sum(x$row.margin)))
  colnames(xmat)<-c(colnames(x$counts),"Total")
 }
 if(!is.na(pct)) {
  xtotal<-sum(x$row.margin)
  dimx<-dim(xmat)
  if(pct[1] == "row") {
   xpctmat<-
    matrix(paste(round(100*x$counts/x$row.margin,ndec),"%",
    sep=""),nrow=dimx[1]-1)
   xpctmat<-rbind(xpctmat,
    paste(round(100*x$col.margin/xtotal,ndec),"%",sep=""))
   xpctmat<-cbind(xpctmat,rep("100%",dimx[1]))
  }
  if(pct[1] == "column") {
   xpctmat<-
    matrix(paste(round(t(t(100*x$counts)/x$col.margin),ndec),"%",
    sep=""),ncol=dimx[2]-1)
   xpctmat<-cbind(xpctmat,
    paste(round(100*x$row.margin/xtotal,ndec),"%",sep=""))
   xpctmat<-rbind(xpctmat,rep("100%",dimx[2]))
  }
  if(pct[1] == "cell") {
   xpctmat<-
    matrix(paste(round(100*x$counts/xtotal,ndec),"%",sep=""),ncol=dimx[2]-1)
   xpctmat<-rbind(xpctmat,
    paste(round(100*x$col.margin/xtotal,ndec),"%",sep=""))
   xpctmat<-cbind(xpctmat,
    paste(c(round(100*x$row.margin/xtotal,ndec),100),"%",sep=""))
  }
  colnames(xpctmat)<-colnames(xmat)
  rownames(xpctmat)<-rownames(xmat)
  if(interdigitate) {
   xpmat<-cbind(xmat[,1],xpctmat[,1])
   for(xcol in 2:dimx[2]) xpmat<-cbind(xpmat,xmat[,xcol],xpctmat[,xcol])
   xpctmat<-xpmat
   colnames(xpctmat)<-
    as.vector(matrix(c(colnames(xmat),rep("%",dimx[2])),
    nrow=2,byrow=TRUE))
  }
 }
 delim.table(xmat,delim=delim,leading.delim=FALSE,label=label)
 if(!is.na(pct)) delim.table(xpctmat,delim=delim,leading.delim=FALSE,label=label)
}
