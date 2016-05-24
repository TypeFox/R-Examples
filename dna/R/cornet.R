cornet = function(data,rescale.scores=FALSE){
 data=as.matrix(data)
 n=as.integer(nrow(data))
 p=as.integer(ncol(data))
 pp=as.integer(p*p) 
 gene.names=colnames(data)
 if (is.null(gene.names))
  gene.names=paste("Gene",1:p)
 out=.C("rcor",as.double(data),s=double(pp),n,p,PACKAGE="dna")
 s=matrix(out$s,p,p,byrow=FALSE) 
 rownames(s)=gene.names
 colnames(s)=gene.names
 if (rescale.scores==TRUE){
  for (i in 1:p) 
   s[i,i]=0
  ss=max(abs(s))
  diag(s)=rep(ss,length(diag(s)))
  s=s/ss
  for (i in 1:p)
   s[i,i]=1
 }
 s
}

