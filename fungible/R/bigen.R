##########################################################
# This function generates binary data with user-defined  #
# thresholds                                             #
#                                                        #
# Arguments:                                             #
# data - population matrix of binary data                #
# n    - number of rows of the generated data            #
# thresholds                                             #
#                                                        #
# Output                                                 #
##########################################################



bigen<-function(data,n,thresholds=NULL, seed=NULL){

if(!is.null(seed)) set.seed(seed)
  
nr<-nrow(data)
nitems<-ncol(data)


##-----------------
## Function rmvnorm by F. Leisch
## a random number generator for the multivariate normal 
## distribution with mean equal to mean and covariance matrix sigma.
rmvnorm<- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean))) 
{
    if (nrow(sigma) != ncol(sigma)) {
        stop("sigma must be a square matrix")
    }
    if (length(mean) != nrow(sigma)) {
        stop("mean and sigma have non-conforming size")
    }
    sigsvd <- svd(sigma)
    retval <- t(sigsvd$v %*% (t(sigsvd$u) * sqrt(sigsvd$d)))
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% retval
    retval <- sweep(retval, 2, mean, "+")
    retval
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## If data supplied, compute population tetrachoric correlation matrix 
if(nr > nitems){
  ## compute thresholds from supplied data
     thresholds<-qnorm(apply(data,2,mean))
     bidat<-matrix(0,nrow=n,ncol=nitems)
     r<-tetcor(X=data)$r
  ## Generate MVN data     
     ghost<-rmvnorm(n, rep(0,nitems), r)
  ## dichotomize data at thresholds   
     for(i in 1:nitems){
         bidat[ghost[,i]<=thresholds[i],i]<-1
     }
    
    result<-list(data=bidat,r=r) 
}

## If population correlation matrix supplied
if(nr==nitems){

     if(is.null(thresholds))stop("thresholds must be supplied with r matrix input")

     bidat<-matrix(0,nrow=n,ncol=nitems)

     if(sum(data) == nitems){
       ghost <- matrix(rnorm(n*nitems),n,nitems)
     }  
     else {  
       ghost<-rmvnorm(n, rep(0,nitems), data)
     }
  ## dichotomize data at thresholds   
     for(i in 1:nitems){
         bidat[ghost[,i]<=thresholds[i],i]<-1
     }   
    result<-list(data=bidat,r=data) 
}

result

}
