intermediate.corr.BO <-
function(n.B, n.O, prop.vec=NULL, prop.list=NULL, corr.vec = NULL, corr.mat=NULL) { 

  validation.bin(n.B, prop.vec)
  validation.ord(n.O, prop.list)

  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }

  n.P<-ncol(corr.mat)-(n.B+n.O)

  cor.mat.BO<-as.matrix(corr.mat[(n.P+1):(n.P+n.B+n.O),(n.P+1):(n.P+n.B+n.O)],(n.B+n.O),(n.B+n.O))

  correlation.bound.check(n.P=0, n.B, n.O, lambda.vec=NULL, prop.vec, prop.list, corr.vec = NULL, corr.mat=cor.mat.BO)
  
  if(n.B==1 && n.O==0) {
  int.mat<-diag(1)
  } else
  if(n.B==0 && n.O==1) {
  int.mat<-diag(1)
  } else
  q.vec=(1-prop.vec)
  marginal<-c(q.vec,prop.list)
  int.mat=ordcont(marginal, Sigma=cor.mat.BO,Spearman = FALSE)$SigmaC

return(int.mat)
}
