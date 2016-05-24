"two2one" <-
function(M,clu=NULL){
	n1<-dim(M)[1]
	n2<-dim(M)[2]
	n<-n1+n2
	M1<-matrix(0,nrow=n,ncol=n)
	M1[1:n1,(n1+1):n]<-M
	dimnames(M1)<-list(unlist(dimnames(M)),unlist(dimnames(M)))
  if(!is.null(clu)) {
    clu<-lapply(clu,function(x)as.numeric(as.factor(x)))
    clu[[2]]<-clu[[2]]+max(clu[[1]])
    clu<-unlist(clu)
  }
 return(list(M=M1,clu=clu))
}

