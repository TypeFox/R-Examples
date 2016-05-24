Int.Corr.NN <-
function(n.NN, corr.vec = NULL, corr.mat=NULL, coef.mat) { 
  
  if(is.null(corr.mat) && !is.null(corr.vec)) {
  d=ceiling(uniroot(function(d) d^2-d-2*length(corr.vec), interval=c(0,1000))$root)
  corr.mat=diag(1,d)
  corr.mat[lower.tri(corr.mat)]=corr.vec
  corr.mat=corr.mat+t(corr.mat)-diag(1,d)
  }

  n.BB<-ncol(corr.mat)-n.NN

  cor.mat.NN<-as.matrix(corr.mat[(n.BB+1):(n.BB+n.NN),(n.BB+1):(n.BB+n.NN)],n.NN,n.NN)

  correlation.bound.check(n.BB=0, n.NN, prop.vec=NULL, corr.vec=NULL, corr.mat=cor.mat.NN, coef.mat)

  if(n.NN==1) {intcor.mat=1
  } else 
  if(n.NN>1){
  usigma.star.c<-diag(n.NN)
  for ( ii in 2:n.NN)        {
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
