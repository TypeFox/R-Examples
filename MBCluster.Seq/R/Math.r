sumRow=function(x,by=NULL){
  if(is.vector(x)) x=matrix(x,ncol=length(x))
  if(is.null(by)) by=rep(1,ncol(x))
  by=as.ordered(as.factor(by))
  mt=outer(by,unique(by),"==")
  xm=x%*%mt
  if(nrow(xm)==1 | ncol(xm)==1) xm=c(xm)
  return(xm)
}

meanRow=function(x,by=NULL){
  if(is.null(by)) by=rep(1,ncol(x))
  by=as.ordered(as.factor(by))
  mt=outer(by,unique(by),"==")
  xm=sweep(x%*%mt,2,sumRow(t(mt)),"/")
  if(nrow(xm)==1 | ncol(xm)==1) xm=c(xm)
  return(xm)
}

sumCol=function(x,by=NULL){
  return(sumRow(t(x),by))
}
meanCol=function(x,by=NULL){
  return(meanRow(t(x),by))
}




maxRow=function(x){
  if(is.vector(x)) return(max(x))
   m=x[,1]
  for(i in 2:ncol(x))
    m=m+(x[,i]-m)*(m<x[,i])
  return(m)
}

minRow=function(x){
  return(-maxRow(-x))
}
