#----------------------------------------------------------------------------------#
# Package: camel                                                                   #
# camel.cmr.mfista(): Multivariate Regression with Calibration                     #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Seq 6th, 2013                                                              #
# Version: 0.2.0                                                                   #
#----------------------------------------------------------------------------------#

camel.cmr.mfista <- function(Y, X, lambda, nlambda, n, d, m, mu, max.ite, prec)
{
  XX = t(X)%*%X
  L = eigen(XX)$values[1]
  beta = array(0,dim=c(nlambda,d,m))
  ite.ext.init = rep(0,nlambda)
  ite.ext.ex = rep(0,nlambda)
  ite.ext.in = rep(0,nlambda)
  str=.C("cmr_sparse_mfista", as.double(Y), as.double(X), 
         as.double(beta), as.integer(n), as.integer(m), as.integer(d), 
         as.double(mu), as.integer(ite.ext.init), as.integer(ite.ext.ex), 
         as.integer(ite.ext.in), as.double(lambda), as.integer(nlambda),
         as.integer(max.ite), as.double(prec), as.double(L), 
         PACKAGE="camel")
  
  beta.list = vector("list", nlambda)
  ite.ext.init = matrix(unlist(str[8]), byrow = FALSE, ncol = nlambda)
  ite.ext.ex = matrix(unlist(str[9]), byrow = FALSE, ncol = nlambda)
  ite.ext.in = matrix(unlist(str[10]), byrow = FALSE, ncol = nlambda)
  ite.ext = vector("list", 3)
  ite.ext[[1]] = ite.ext.init
  ite.ext[[2]] = ite.ext.ex
  ite.ext[[3]] = ite.ext.in
  dm = d*m
  for(i in 1:nlambda){
    beta.i = matrix(unlist(str[3])[((i-1)*dm+1):(i*dm)],nrow=d,ncol=m,byrow = FALSE)
    beta.list[[i]] = beta.i
  }
  
  return(list(beta=beta.list, ite=ite.ext))
}
