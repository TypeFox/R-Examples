intermediate.corr.PC <-
function(n.P, n.B, n.C, lambda.vec=NULL, prop.vec=NULL, coef.mat=NULL, corr.vec = NULL, corr.mat=NULL){

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
 
  if(n.C>0) {
  samples = 1e+05
  xmat=matrix(NA, nrow=samples, ncol=n.C)
  for (i in 1:n.C){
  #x=as.vector(rnorm(samples))
  x=as.vector(n1)
  xx=cbind(1,x,x^2,x^3)
  xmat[,i]=xx%*%coef.mat[,i]
  }
  mydata=xmat
  } 
    
  amat=matrix(rep(as.vector(cor(mydata,n2)),n.P),n.P,n.C, byrow=T)
  
  corr.mat.PN=matrix(c(corr.mat[1:n.P,(n.P+n.B+1):(n.P+n.B+n.C)]),n.P,n.C)
 
  corr.mat.PN2=matrix(0,n.P,n.C)
  for(i in 1:n.P) {
  for(j in 1:n.C) {
  corr.mat.PN2[i,j]=corr.mat.PN[i,j]/amat[i,j]
  }
  }

  u=runif(samples)
  zu=matrix(qnorm(u),samples,1)
  pu=sapply(1:n.P, function(i) qpois(u,lambda.vec[i]))

  chat=matrix(rep(cor(pu,zu),n.C),n.P,n.C)

  int.corrmat=matrix(0,n.P,n.C)
  for(i in 1:n.P){
  for(j in 1:n.C) {
  int.corrmat[i,j]=corr.mat.PN2[i,j]/chat[i,j]
  }
  }

return(int.corrmat)
}
