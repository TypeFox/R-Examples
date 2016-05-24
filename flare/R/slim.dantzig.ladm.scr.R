#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# slim.dantzig.ladm.scr(): Regression with Dantzig()                               #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jul 8th, 2014                                                              #
# Version: 1.4.0                                                                   #
#----------------------------------------------------------------------------------#

slim.dantzig.ladm.scr <- function(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
{
  if(verbose==TRUE)
    cat("Dantzig selector with screening.\n")
  XY = crossprod(X,Y/n)
  XX = crossprod(X,X/n)#+0.05*lambda[nlambda]*diag(d)
  beta = matrix(0,nrow=d,ncol=nlambda)
  ite.int = rep(0,nlambda)
  ite.int1 = rep(0,nlambda)
  ite.int2 = rep(0,nlambda)
  if(intercept) {
    intcep=1
  }else{
    intcep=0
  }
  
  if(n<=3){
    num.scr1 = n
    num.scr2 = n
  }else{
    num.scr1 = ceiling(n/log(n))
    num.scr2 = n-1
  }

  order0 = order(abs(XY),decreasing = TRUE)
  idx.scr = order0; num.scr = length(idx.scr)
  idx.scr1 = order0[1:num.scr1]
  idx.scr2 = order0[1:num.scr2]
  X1 = X[,idx.scr]
  XXX = crossprod(X1,crossprod(tcrossprod(X1,X1),X1))/(n^2)
  gamma = max(colSums(abs(XXX)))
  str=.C("slim_dantzig_ladm_scr", as.double(XY), as.double(XX), as.double(XXX), 
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
