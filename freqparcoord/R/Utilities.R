
# author:  Norm Matloff

#####################################################################
# add positive jitter
#####################################################################

# similar to jitter(), but only generating postive values in (0,1);
# typical example is age, which in many data sets is truncated to the
# lowest integer

posjitter <- function(x) x + runif(length(x))

#####################################################################
# generate random variates from mixtures of MV normals
#####################################################################

# arguments:
#    n:  number of random vectors to generate
#    dm:  length of each vector
#    nmix:  number of MV normal distributions in the mixturee
#    means:  mean vectors of the MV normal distributions; an R list of nmix
#            vectors of length dm each
#    covs:  covariance matrices of the MV normal distributions; an R list of nmi
#           matrices, each dm x dm 
#    wts:  mixture probabilities

# return value: 

#   n random vectors of length dm, stored in an n x dm matrix; grouped
#   by MV normal distribution 

# example:

#    cv <- 0.5*diag(2)
#    x <- rmixmvnorm(25000,2,3,list(c(0,0),c(1,2),c(3,3)),list(cv,cv,cv))
#    smoothScatter(x)

rmixmvnorm <- function(n,dm,nmix,means,covs,wts=rep(1/nmix,nmix)) {
   labels <- sample(1:nmix,n,replace=T,prob=wts) 
   lt <- table(labels)
   avvlbls <- as.integer(names(lt))
   out <- matrix(nrow=n,ncol=dm)
   rownum <- 1
   for (i in avvlbls) {
      ni <- lt[i]
      rows <- rownum:(rownum+ni-1)
      out[rows,] <- rmvnorm(ni,mean=means[[i]],sigma=covs[[i]])
      rownum <- rownum + ni
   }
   out
}
