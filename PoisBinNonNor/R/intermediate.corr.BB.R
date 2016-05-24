intermediate.corr.BB <-
function(n.P,n.B,n.C, prop.vec, corr.vec = NULL, corr.mat=NULL) { 

  validation.bin(n.B, prop.vec)

  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }

  if( (n.P+n.B+n.C) !=ncol(corr.mat)) stop("Dimension of the correlation matrix is misspecied!")

  cor.mat.BB<-as.matrix(corr.mat[(n.P+1):(n.P+n.B),(n.P+1):(n.P+n.B)],n.B,n.B)

  correlation.bound.check(n.P=0, n.B, n.C=0, lambda.vec=NULL, prop.vec, coef.mat=NULL, corr.vec = NULL, corr.mat=cor.mat.BB)

  q.vec=(1-prop.vec)

  if(n.B==1) {tetcor.mat=1
  } else 
  if(n.B>1) {
  ##get the pairwise correlations##
  usigma.star.b<-diag(n.B)
  for ( i  in 2:n.B)            {
  for ( ii in 1:(i-1))          {
  mycorfuncb<- function(ro.mat) {
  r<-integrate(function(z2) {sapply(z2, function(z2) {integrate(function(z1) ((2*pi*sqrt((1-ro.mat^2)))^-1)* exp(-(z1^2-2*ro.mat*z1*z2+z2^2)/(2*(1-ro.mat^2))), -Inf, qnorm(prop.vec[i]) )$value})},  -Inf, qnorm(prop.vec[ii]))$value-
  cor.mat.BB[ii,i]*sqrt(prop.vec[i]*q.vec[ii]*prop.vec[ii]*q.vec[i])-(prop.vec[i]*prop.vec[ii])
  r
  }
  p0 <-0
  usigma.star.b[ii,i]=suppressWarnings(dfsane(par = p0, fn=mycorfuncb, control=list(trace=FALSE)))$par
  }#ii
  }#i
  tetcor.mat<-usigma.star.b+t(usigma.star.b)
  diag(tetcor.mat)<-1
  }#if

return(tetcor.mat)
}
