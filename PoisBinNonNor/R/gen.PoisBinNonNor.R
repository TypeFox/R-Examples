gen.PoisBinNonNor<-function(n, n.P, n.B, n.C, lambda.vec=NULL, prop.vec=NULL, mean.vec=NULL, variance.vec=NULL, coef.mat=NULL, final.corr.mat){

   if(missing(n)==TRUE)          stop("n was not specified! \n")
   if(missing(final.corr.mat))   stop("Final correlation matrix was not specified! \n")

   validation.bin(n.B, prop.vec)

   if(ncol(final.corr.mat) != (n.P+n.B+ n.C)) {
     stop("Dimension of final correlation matrix does not match the number of variables! \n")
   } #if

   myz<-rmvnorm(n, mean=rep(0,(n.P+n.B+n.C)),final.corr.mat)
   samples=1e+05

   if(!is.null(lambda.vec) && is.null(prop.vec) && is.null(coef.mat) )  {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]))),n,n.P)
   mydata=cbind(myb1)
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && is.null(coef.mat) )  {
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   mydata=cbind(myb2)
   } else
   if(is.null(lambda.vec) && is.null(prop.vec) && !is.null(coef.mat) )  {
   myb3<-matrix(0,n,n.C)
   myb3=matrix(sapply(1:n.C, function(j) sapply(1:n, function(i)   (coef.mat[1,j]+coef.mat[2,j]*myz[i,j]+coef.mat[3,j]*(myz[i,j]^2)+coef.mat[4,j]*myz[i,j]^3)*
   sqrt(variance.vec[j])+(mean.vec[j]))),n,n.C)
   mydata=cbind(myb3)
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && is.null(coef.mat) ) {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]) )),n,n.P)
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii+n.P]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   mydata=cbind(myb1,myb2)
   } else
   if(!is.null(lambda.vec) && is.null(prop.vec) && !is.null(coef.mat) ) {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]))),n,n.P)
   myb3<-matrix(0,n,n.C)
   myb3=matrix(sapply(1:n.C, function(j) sapply(1:n, function(i)   (coef.mat[1,j]+coef.mat[2,j]*myz[i,j+n.P]+coef.mat[3,j]*(myz[i,j+n.P]^2)+coef.mat[4,j]*myz[i,j+n.P]^3)*
   sqrt(variance.vec[j])+(mean.vec[j]))),n,n.C)
   mydata=cbind(myb1,myb3)
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && !is.null(coef.mat) ) {
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii+n.P]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   myb3<-matrix(0,n,n.C)
   myb3=matrix(sapply(1:n.C, function(j) sapply(1:n, function(i)   (coef.mat[1,j]+coef.mat[2,j]*myz[i,j+n.B]+coef.mat[3,j]*(myz[i,j+n.B]^2)+coef.mat[4,j]*myz[i,j+n.B]^3)*
   sqrt(variance.vec[j])+(mean.vec[j]))),n,n.C)
   mydata=cbind(myb2,myb3)
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && !is.null(coef.mat) ) {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]) )),n,n.P)
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii+n.P]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   myb3<-matrix(0,n,n.C)
   myb3=matrix(sapply(1:n.C, function(j) sapply(1:n, function(i)   (coef.mat[1,j]+coef.mat[2,j]*myz[i,j+n.P+n.B]+coef.mat[3,j]*(myz[i,j+n.P+n.B]^2)+coef.mat[4,j]*myz[i,j+n.P+n.B]^3)*
   sqrt(variance.vec[j])+(mean.vec[j]))),n,n.C)
   mydata=cbind(myb1,myb2,myb3)
   }
colnames(mydata)=NULL
return(mydata)
}