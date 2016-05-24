#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# slim.dantzig.ladm.scr2(): Regression with Dantzig()                              #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jul 8th, 2014                                                              #
# Version: 1.4.0                                                                   #
#----------------------------------------------------------------------------------#

slim.dantzig.ladm.scr2 <- function(Y, X, lambda, nlambda, n, d, maxdf, rho, max.ite, prec, intercept, verbose)
{
  if(verbose==TRUE)
    cat("Dantzig selector with screening.\n")
  XY = crossprod(X,Y)/n
  XX = crossprod(X)/n#+0.05*lambda[nlambda]*diag(d)
  beta = matrix(0,nrow=d,ncol=nlambda)
  ite.int = rep(0,nlambda)
  ite.int1 = rep(0,nlambda)
  if(intercept) {
    intcep=1
  }else{
    intcep=0
  }
  if(d<=3){
    num.scr1 = d
  }else{
    num.scr1 = ceiling(d/log(d))
  }
  order0 = order(abs(XY),decreasing = TRUE)
  idx.scr = order0; num.scr = length(idx.scr)
  idx.scr1 = order0[1:num.scr1]
  XX1 = XX[idx.scr,idx.scr]
  XXX = crossprod(XX1,XX1)
  gamma = max(colSums(abs(X)))
  begt = Sys.time()
  str=.C("slim_dantzig_ladm_scr2", as.double(XY), as.double(XX), as.double(XXX), 
         as.double(beta), as.integer(n), as.integer(d), as.double(rho),
         as.integer(ite.int), as.integer(ite.int1), as.integer(num.scr1), 
         as.integer(idx.scr), as.integer(idx.scr1), as.double(gamma), as.double(lambda), 
         as.integer(nlambda), as.integer(max.ite), as.double(prec), 
         as.integer(intcep),PACKAGE="flare")
  runt = Sys.time() - begt
  beta.list = vector("list", nlambda)
  for(i in 1:nlambda){
    beta.i = unlist(str[4])[((i-1)*d+1):(i*d)]
    beta.list[[i]] = beta.i
  }
  ite.int = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  ite.int1 = matrix(unlist(str[9]), byrow = FALSE, ncol = nlambda)
  ite = list()
  ite[[1]] = ite.int1
  ite[[2]] = ite.int
  return(list(beta=beta.list, ite=ite))
}
