PCnet = function(data, ncom=3, rescale.data=TRUE, symmetrize.scores=TRUE,
rescale.scores=FALSE){
 data=as.matrix(data)
 n=as.integer(nrow(data))
 p=as.integer(ncol(data))
 pp=as.integer(p*p) 
 gene.names=colnames(data)
 if (is.null(gene.names))
  gene.names=paste("Gene",1:p)
 out=.C("rpcnet",as.double(data),s=double(pp),as.integer(ncom),n,p,
as.integer(rescale.data),as.integer(symmetrize.scores),
as.integer(rescale.scores),PACKAGE="dna")
 s=matrix(out$s,p,p,byrow=FALSE) 
 rownames(s)=gene.names
 colnames(s)=gene.names
 s
}

