###########################################################################################################################
###Simulates a sample of size n from a set of multivariate Poisson, binary and ordinal variables.
###########################################################################################################################

gen.PoisBinOrd<-function(n, n.P, n.B, n.O, lambda.vec=NULL, prop.vec=NULL, prop.list=NULL, final.corr.mat){

   if(missing(n)==TRUE)          stop("n was not specified! \n")
   if(missing(final.corr.mat))   stop("Final correlation matrix was not specified! \n")

   validation.bin(n.B, prop.vec)
   validation.ord(n.O, prop.list)


   if(ncol(final.corr.mat) != (n.P+n.B+n.O)) {
     stop("Dimension of final correlation matrix does not match the number of variables! \n")
   } #if

   myz<-rmvnorm(n, mean=rep(0,(n.P+n.B+n.O)),final.corr.mat)
   samples=1e+05

   if(!is.null(lambda.vec) && is.null(prop.vec) && is.null(prop.list) )  {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]) )),n,n.P)
   mydata=cbind(myb1)
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && is.null(prop.list) )  {
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   mydata=cbind(myb2)
   } else
   if(is.null(lambda.vec) && is.null(prop.vec) && !is.null(prop.list) )  {
   #myb3=ordsample(n, marginal=prop.list, Sigma=final.corr.mat, Spearman=FALSE, cormat="continuous")
   myb3<-matrix(0,samples,n.O)
   for(k in 1:n.O) {
   cv=qnorm(prop.list[[k]])
   for(i in 1:samples){
   myb3[i,k]=length(which(cv<myz[i,k]))+1
   }
   } 
   mydata=cbind(myb3)
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && is.null(prop.list) ) {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]) )),n,n.P)
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii+n.P]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   mydata=cbind(myb1,myb2)
   } else
   if(!is.null(lambda.vec) && is.null(prop.vec) && !is.null(prop.list) ) {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]))),n,n.P)
   #myb3=ordsample(n, marginal=prop.list, Sigma=final.corr.mat[(n.P+1):(n.P+n.O),(n.P+1):(n.P+n.O)], Spearman=FALSE, cormat="continuous")
   myb3<-matrix(0,samples,n.O)
   for(k in 1:n.O) {
   cv=qnorm(prop.list[[k]])
   for(i in 1:samples){
   myb3[i,k]=length(which(cv<myz[i,k+(n.P+n.B)]))+1
   }
   }
   mydata=cbind(myb1,myb3)
   } else
   if(is.null(lambda.vec) && !is.null(prop.vec) && !is.null(prop.list) ) {
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii+n.P]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   #myb3=ordsample(n, marginal=prop.list, Sigma=final.corr.mat[(n.P+1):(n.P+n.O),(n.P+1):(n.P+n.O)], Spearman=FALSE, cormat="continuous")
   myb3<-matrix(0,samples,n.O)
   for(k in 1:n.O) {
   cv=qnorm(prop.list[[k]])
   for(i in 1:samples){
   myb3[i,k]=length(which(cv<myz[i,k+(n.P+n.B)]))+1
   }
   }
   mydata=cbind(myb2,myb3)
   } else
   if(!is.null(lambda.vec) && !is.null(prop.vec) && !is.null(prop.list) ) {
   myb1=matrix(sapply(1:n.P, function(ii) sapply(1:n, function(i)  qpois(pnorm(myz[i,ii]),lambda.vec[ii]) )),n,n.P)
   myb2<-matrix(0,n,n.B)
   myb2=matrix(sapply(1:n.B, function(ii) sapply(1:n, function(i)  if(1*myz[i,ii+n.P]>qnorm(1-prop.vec[ii])) myb2[i,ii]=1 else myb2[i,ii]=0 )),n,n.B)
   #myb3=ordsample(n, marginal=prop.list, Sigma=final.corr.mat[(n.P+n.B+1):(n.P+n.B+n.O),(n.P+n.B+1):(n.P+n.B+n.O)], Spearman=FALSE, cormat="continuous")
   myb3<-matrix(0,samples,n.O)
   for(k in 1:n.O) {
   cv=qnorm(prop.list[[k]])
   for(i in 1:samples){
   myb3[i,k]=length(which(cv<myz[i,k+(n.P+n.B)]))+1
   }
   }
   mydata=cbind(myb1,myb2,myb3)
   }

return(mydata)
}
