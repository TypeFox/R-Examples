#-------------------------------------------------------------------------------#
# Package: fastclime                                                            #
# fastclime(): Main Function                                                    #
# Authors: Haotian Pang, Di Qi, Han Liu and Robert Vanderbei                    #
# Emails: <hpang@princeton.edu>, <dqi@princeton,edu>, <hanliu@princeton.edu>    #
# and <rvdb@princetonedu>                                                       #
# Date: April 22th 2016                                                         #
# Version: 1.4.1          					                                            #
#-------------------------------------------------------------------------------#
fastclime <- function(x, lambda.min=0.1, nlambda = 50)
{
  
  gcinfo(FALSE)
  cov.input<-1
  SigmaInput<-x
  d<-dim(SigmaInput)[2]

  if(!isSymmetric(x))
  {
     n<-dim(SigmaInput)[1]
     SigmaInput<-cov(x)*(1-1/n)
     cov.input<-0
  }


 
  cat("Allocating memory \n")
  maxnlambda=0
  mu_input<-matrix(0,nlambda,d)
  iicov<-matrix(0,nlambda,d*d)
  lambdamin<-lambda.min

 cat("start recovering \n")  
     str=.C("parametric", as.double(SigmaInput), as.integer(d), as.double(mu_input), as.double(lambdamin), as.integer(nlambda), as.integer(maxnlambda), as.double(iicov), PACKAGE="fastclime")

 cat("preparing precision and path matrix list \n")
 
  sigmahat<-matrix(unlist(str[1]),d)   
  mu<-matrix(unlist(str[3]),nlambda,d)
  maxnlambda<-unlist(str[6])+1
  iicov<-matrix(unlist(str[7]),nlambda,d*d)
  mu<-mu[1:maxnlambda,]
  icov<-list()
  
  for (i in 1:maxnlambda)
  {
     icov[i]<-list(matrix(iicov[i,],d,d))
  }
  #icov[maxnlambda+1]=list(icov[[maxnlambda]])

  result<-list("data" = x, "cov.input" = cov.input, "sigmahat" = sigmahat, "maxnlambda" = maxnlambda, "lambdamtx" = mu, "icovlist" = icov)

  rm(x,cov.input,sigmahat,maxnlambda,mu,icov,iicov,
    nlambda,lambdamin,mu_input,SigmaInput,d)
  gc()
  class(result) = "fastclime"
  cat("Done! \n")
  return(result)

}

print.fastclime = function(x, ...)
{	

	if(x$cov.input) cat("Input: The Covariance Matrix\n")
	if(!x$cov.input) cat("Input: The Data Matrix\n")
	cat("Path length:",x$nlambda,"\n")
	cat("Graph dimension:",ncol(x$data),"\n")
	#cat("Sparsity level:",min(x$sparsity),"----->",max(x$sparsity),"\n")
}


plot.fastclime = function(x, ...){
        gcinfo(FALSE)
        s<-x$lambda[,1]
        poslambda<-s[s>0]
        npos<-length(poslambda)
	
	plot(x$lambda[1:npos,1], x$sparsity[1:npos], log = "x", xlab = "Regularization Parameter", ylab = "Sparsity Level", type = "l", main = "Sparsity vs. Regularization")
	

}
