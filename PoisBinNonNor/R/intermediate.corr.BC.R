intermediate.corr.BC <-
function(n.P, n.B, n.C, lambda.vec=NULL, prop.vec=NULL, coef.mat=NULL, corr.vec = NULL, corr.mat=NULL) { 

  validation.bin(n.B, prop.vec)

  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }

  if((n.P+n.B+n.C) !=ncol(corr.mat)) stop("Dimension of the correlation matrix is misspecied!")

  correlation.bound.check(n.P, n.B, n.C, lambda.vec, prop.vec, coef.mat, corr.vec = NULL, corr.mat=corr.mat)

  cor.bn=apply(as.matrix(prop.vec,length(prop.vec),1),1,function(x) dnorm(qnorm(x))/sqrt(x*(1-x)))
  bicor.mat=t(sapply(1:n.B, function(ii) sapply(1:n.C, function(i) (corr.mat[(ii+n.P),(i+n.B+n.P)]/cor.bn[ii])/(coef.mat[2,i]+3*coef.mat[4,i]))))

return(bicor.mat)
}
