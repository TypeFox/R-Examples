###########################################################################################################################
###Simulates a sample of size n from a set of multivariate binary and nonnormal continuous variables.
###########################################################################################################################

gen.Bin.NonNor<-function(n, n.BB, n.NN, prop.vec=NULL, mean.vec=NULL, variance.vec=NULL, skewness.vec=NULL, kurtosis.vec=NULL, final.corr.mat, coef.mat=NULL){

  if(missing(n)==TRUE)          stop("n was not specified! \n")
  if(missing(final.corr.mat))   stop("Final correlation matrix was not specified! \n")

  validation.bin(n.BB, prop.vec)
  validation.skewness.kurtosis(n.NN, skewness.vec, kurtosis.vec)

  if(n.BB >0 && n.NN >0 && ncol(final.corr.mat) != length(prop.vec)+ length(skewness.vec)) {
     stop("Dimension of final correlation matrix does not match the number of variables! \n")
  }  else
  if(n.BB ==0 && n.NN >0 && ncol(final.corr.mat) != length(skewness.vec)) {
     stop("Dimension of final correlation matrix does not match the number of continuous variables! \n")
  }  else
  if(n.BB >0 && n.NN==0 && ncol(final.corr.mat) != length(prop.vec)) {
     stop("Dimension of final correlation matrix does not match the number of binary variables! \n")
  }#if

  if(!is.null(prop.vec)&& is.null(skewness.vec)) {
  myz<-rmvnorm(n, mean=rep(0,n.BB),final.corr.mat)
  myb<-matrix(0,n,n.BB)
  myb=matrix(sapply(1:n.BB, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii]>qnorm(1-prop.vec[ii])) myb[i,ii]=1 else myb[i,ii]=0 )),n,n.BB)
  mydata=cbind(myb)
  } else
  if(is.null(prop.vec) && !is.null(skewness.vec)) {
  myz<-rmvnorm(n, mean=rep(0,n.NN),final.corr.mat)
  myy<-matrix(0,n,n.NN)
  myy=matrix(sapply(1:n.NN, function(j) sapply(1:n,function(i)   (coef.mat[1,j]+coef.mat[2,j]*myz[i,j]+coef.mat[3,j]*(myz[i,j]^2)+coef.mat[4,j]*myz[i,j]^3)*
  sqrt(variance.vec[j])+(mean.vec[j]))),n,n.NN)
  mydata=cbind(myy)
  }else
  if(!is.null(prop.vec) && !is.null(skewness.vec))  {
  myz<-rmvnorm(n, mean=rep(0,(n.BB+n.NN)),final.corr.mat)
  myb<-matrix(0,n,n.BB)
  myb=matrix(sapply(1:n.BB, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii]>qnorm(1-prop.vec[ii])) myb[i,ii]=1 else myb[i,ii]=0)),n,n.BB)
  myb=myb
  myy<-matrix(0,n,n.NN)
  myy=matrix(sapply(1:n.NN, function(j) sapply(1:n,function(i)   (coef.mat[1,j]+coef.mat[2,j]*myz[i,j+n.BB]+coef.mat[3,j]*(myz[i,j+n.BB]^2)+coef.mat[4,j]*myz[i,j+n.BB]^3)*
  sqrt(variance.vec[j])+(mean.vec[j]))),n,n.NN)
  mydata=cbind(myb,myy)
  colnames(mydata)<-NULL
  }#if

return(mydata)
}
