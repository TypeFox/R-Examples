intermediate.corr.PP <-
function(n.P, lambda.vec, corr.vec=NULL, corr.mat=NULL){

  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }

  cor.mat.P<-as.matrix(corr.mat[1:n.P,1:n.P],n.P,n.P)

  correlation.bound.check(n.P, n.B=0, n.O=0, lambda.vec, prop.vec=NULL, prop.list=NULL, corr.vec = NULL, corr.mat=cor.mat.P)

  samples=1e+05
  if(n.P==1) {
  int.mat<-diag(1)
  } else 
  if(n.P>1) {
  umat=runif(samples)
  xmat1=sapply(1:length(lambda.vec), function(i) qpois(umat,lambda.vec[i]))
  upp.lim=cor(xmat1)[col(cor(xmat1)) > row(cor(xmat1))] 
  xmat2=sapply(1:length(lambda.vec),function(i) qpois(1-umat,lambda.vec[i]))
  low.lim=cor(xmat1,xmat2)[col(cor(xmat1,xmat2)) > row(cor(xmat1,xmat2))]   
    
  anymat=diag(1,n.P)
  anymat[upper.tri(anymat)]=upp.lim
  anymat[lower.tri(anymat)]=low.lim

  avec=-(low.lim*upp.lim)/(low.lim+upp.lim)
  amat=diag(1,n.P)
  amat[upper.tri(amat)]=avec
  #amat[lower.tri(amat)]=avec
  amat=amat+t(amat)
  diag(amat)=1

  int.mat<-matrix(0,n.P,n.P)
  for(i  in 2:n.P){
  for(j  in 1:(i-1)){
  int.mat[j,i]<-log((cor.mat.P[j,i]+amat[j,i])/amat[j,i])/log((anymat[j,i]+amat[j,i])/amat[j,i])
  }
  }
  int.mat<-int.mat+t(int.mat)
  diag(int.mat)=1
  }#if

return(int.mat)
}
