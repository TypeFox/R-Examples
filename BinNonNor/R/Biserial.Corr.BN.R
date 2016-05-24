Biserial.Corr.BN <-
function(n.BB, n.NN, prop.vec, corr.vec = NULL, corr.mat=NULL, coef.mat) { 

  validation.bin(n.BB, prop.vec)

  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }

  correlation.bound.check(n.BB, n.NN, prop.vec, corr.vec=NULL, corr.mat, coef.mat)

  cor.bn=apply(as.matrix(prop.vec,length(prop.vec),1),1,function(x) dnorm(qnorm(x))/sqrt(x*(1-x)))
  bicor.mat=t(sapply(1:n.BB, function(ii) sapply(1:n.NN, function(i) (corr.mat[ii,i+n.BB]/cor.bn[ii])/(coef.mat[2,i]+3*coef.mat[4,i]))))

return(bicor.mat)
}
