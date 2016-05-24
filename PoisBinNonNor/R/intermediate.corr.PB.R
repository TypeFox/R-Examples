intermediate.corr.PB <-
function(n.P, n.B, n.C, lambda.vec=NULL, prop.vec=NULL, coef.mat=NULL, corr.vec = NULL, corr.mat=NULL){

  validation.bin(n.B, prop.vec)

  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }
 
  if((n.P+n.B+n.C) !=ncol(corr.mat)) stop("Dimension of the correlation matrix is misspecied!")

  correlation.bound.check(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec = NULL, corr.mat=corr.mat)
  
  samples=1e+05
  n1=rnorm(samples)
  n2=n1
 
  if(n.B>0) {
  mydatabin=matrix(0,samples,n.B)
  for(k in 1:n.B){
  cv=qnorm(1-prop.vec[k])
  for(i in 1:samples){
  if(n1[i]>cv) mydatabin[i,k]=1
  }
  }
  mydata=mydatabin
  } 
    
  amat=matrix(rep(as.vector(cor(mydata,n2)),n.P),n.P,n.B, byrow=T)
  corr.mat.PO=matrix(c(corr.mat[1:n.P,(n.P+1):(n.P+n.B)]),n.P,n.B)
  corr.mat.PN2=matrix(0,n.P,n.B)

  for(i in 1:n.P){
  for(j in 1:n.B) {
  corr.mat.PN2[i,j]=corr.mat.PO[i,j]/amat[i,j]
  }
  }

  u=runif(samples)
  zu=matrix(qnorm(u),samples,1)
  pu=sapply(1:n.P, function(i) qpois(u,lambda.vec[i]))

  chat=matrix(rep(cor(pu,zu),n.B),n.P,n.B)

  int.corrmat=matrix(0,n.P,n.B)
  for(i in 1:n.P){
  for(j in 1:n.B) {
  int.corrmat[i,j]=corr.mat.PN2[i,j]/chat[i,j]
  }
  }

return(int.corrmat)
}
