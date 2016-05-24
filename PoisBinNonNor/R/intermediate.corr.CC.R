intermediate.corr.CC <-
function(n.P, n.B, n.C, coef.mat=NULL, corr.vec = NULL, corr.mat=NULL) { 
  
  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }

  if((n.P+n.B+n.C) !=ncol(corr.mat)) stop("Dimension of the correlation matrix is misspecied!")

  cor.mat.NN<-as.matrix(corr.mat[(n.P+n.B+1):(n.P+n.B+n.C),(n.P+n.B+1):(n.P+n.B+n.C)],n.C,n.C)

  correlation.bound.check(n.P=0, n.B=0, n.C, lambda.vec=NULL, prop.vec=NULL, coef.mat, corr.vec = NULL, corr.mat=cor.mat.NN)
 
  if(n.C==1) {intcor.mat=1
  } else 
  if(n.C>1){
  usigma.star.c<-diag(n.C)
  for ( ii in 2:n.C)         {
  for ( i in 1:(ii-1))       {
  mycorfunc<- function(rro)  {
  r <- rro*((coef.mat[2,i]*coef.mat[2,ii])+(3*coef.mat[2,i]*coef.mat[4,ii])+(3*coef.mat[4,i]*coef.mat[2,ii])+(9*coef.mat[4,i]*coef.mat[4,ii]))+
       (rro^2)*(2*coef.mat[3,i]*coef.mat[3,ii])+ rro^3*(6*coef.mat[4,i]*coef.mat[4,ii])-cor.mat.NN[i,ii]
  r
  }#myfunc
  p0 <-0
  usigma.star.c[i,ii]=suppressWarnings(dfsane(par = p0, fn=mycorfunc, control=list(trace=FALSE)))$par
  }#ii
  }#i
  intcor.mat<-usigma.star.c+t(usigma.star.c)
  diag(intcor.mat)<-1
  }#if 

return(intcor.mat)
}
