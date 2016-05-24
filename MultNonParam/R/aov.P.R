aov.P<-function(dattab,permi=NULL,be=NULL){
# Rows represent treatments, columns represent blocks
  if(is.matrix(dattab)){
     n<-prod(dim(dattab))
     nb<-dim(dattab)[2]
     ng<-dim(dattab)[1]
     permi<-rep(1:ng,nb)
     be<-seq(nb)*ng
     xx<-as.vector(dattab)
  }else{
     n<-length(dattab)
     xx<-dattab
     if(is.null(be)) be<-n
     ng<-length(unique(permi))
     nb<-length(be)
  }
  out<-.Fortran("aovp",
                as.integer(n),
                as.integer(permi),
                as.integer(nb),
                as.integer(be),
                as.double(xx),
                tot=as.double(0),pv=as.double(0)
                ,PACKAGE="MultNonParam")
  return(list(pv=out$pv,tot=out$tot))
}
