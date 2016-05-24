
"dim.repweights_compressed"<-function(x){
  c(length(x$index),ncol(x$weights))
}

"dimnames.repweights_compressed"<-function(x){
  list(names(x$index), colnames(x$weights))
}

"[.repweights_compressed"<-function(x,i,...,drop=FALSE){
  if (!missing(i)){
    x$index<-x$index[i]
    if(!missing(..1))
      x$weights<-x$weights[,..1,drop=FALSE]
  } else{
      ## this is faster than just subscripting x$weights (!)
      x<-list(index=x$index,
              weights=x$weights[,...,drop=FALSE])
      class(x)<-c("repweights_compressed","repweights")
  }
  x
}

"as.matrix.repweights_compressed"<-function(x,...){
  x$weights[x$index,,drop=FALSE]
}

"as.vector.repweights_compressed"<-function(x,...){
  as.vector(x$weights[x$index,])
}

"as.matrix.repweights"<-function(x,...){
  x
}

compressWeights<-function(rw,...){
  UseMethod("compressWeights")
}

"compressWeights.repweights_compressed"<-function(rw,...){
  compressWeights(as.matrix(rw))
}

compressWeights.default<-function(rw,...){
  mat<-as.matrix(rw)
  tmp<-apply(mat,1,function(x) paste(x,collapse="\r"))
  unq<-!duplicated(mat)
  rval<-list(weights=mat[unq,],index=match(tmp,tmp[unq]))
  class(rval)<-c("repweights_compressed","repweights")
  rval
}

compressWeights.svyrep.design<-function(rw,...){
  rw$repweights<-compressWeights(rw$repweights,...)
  rw
}
