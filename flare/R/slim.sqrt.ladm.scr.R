#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# slim.sqrt.ladm.scr(): Regression with sqrt lasso()                               #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jul 8th, 2014                                                              #
# Version: 1.4.0                                                                   #
#----------------------------------------------------------------------------------#

slim.sqrt.ladm.scr <- function(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
{
  if(verbose==TRUE)
    cat("Square root Lasso with screening.\n")
  XY = crossprod(X,Y)
  beta = matrix(0,nrow=d,ncol=nlambda)
  ite.int = rep(0,nlambda)
  ite.int1 = rep(0,nlambda)
  ite.int2 = rep(0,nlambda)
  if(intercept) {
    intcep=1
  }else{
    intcep=0
  }
  if(d>n){
    if(n<=3){
      num.scr1 = n
      num.scr2 = n
    }else{
      num.scr1 = ceiling(n/log(n))
      num.scr2 = n-1
    }
  }else{
    if(d<=3){
      num.scr1 = d
      num.scr2 = d
    }else{
      num.scr1 = ceiling(sqrt(d))
      num.scr2 = ceiling(d/log(d))
    }
  }
  order0 = order(XY,decreasing = TRUE)
  idx.scr = order0; num.scr = length(idx.scr)
  idx.scr1 = order0[1:num.scr1]
  idx.scr2 = order0[1:num.scr2]
  X1 = X[,idx.scr]
  XX = crossprod(X1)
  gamma = max(colSums(abs(XX)))
  str=.C("slim_sqrt_ladm_scr", as.double(Y), as.double(X), as.double(XX), 
         as.double(beta), as.integer(n), as.integer(d), as.double(rho),
         as.integer(ite.int), as.integer(ite.int1), as.integer(ite.int2), 
         as.integer(num.scr1), as.integer(num.scr2), 
         as.integer(idx.scr), as.integer(idx.scr1), as.integer(idx.scr2), 
         as.double(gamma), as.double(lambda), 
         as.integer(nlambda), as.integer(max.ite), as.double(prec), 
         as.integer(intcep),PACKAGE="flare")
  beta.list = vector("list", nlambda)
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
  }
  ite.int = unlist(str[8])
  ite.int1 = unlist(str[9])
  ite.int2 = unlist(str[10])
  ite = list()
  ite[[1]] = ite.int1
  ite[[2]] = ite.int2
  ite[[3]] = ite.int
  return(list(beta=beta.list, ite=ite))
}
