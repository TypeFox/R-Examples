linesROC <-
function(trueMatrix,pep,col="red",lty=1,lwd=1) {
  n<-sum(!is.na(trueMatrix))
  l<-sort(pep,decreasing=TRUE)
  pos<-which(trueMatrix==1)
  npos<-length(pos)
  neg<-which(trueMatrix==0)
  nneg<-length(neg)
  A<-NA
  if (nneg==0 || npos==0) {
    cat("ROC curve undefined.\n")
  } else {
    x0<-0
    y0<-0
    A<-0
    for (i in 1:length(l)) {
      x1<-sum(1*(pep[neg]>=l[i]))/nneg
      y1<-sum(1*(pep[pos]>=l[i]))/npos
      if (x1>x0) {A<-A+(x1-x0)*(y0+y1)/2}
      segments(x0,y0,x1,y1,col=col,lty=lty)
      x0<-x1
      y0<-y1
    }
  }
  cat(col,"ROC area =",A,"\n")
  return(A)
}
