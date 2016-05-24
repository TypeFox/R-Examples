#----------------------------------------------------------------------------------#
# Package: camel                                                                   #
# camel.tiger.clime.hadm(): Coordinate descent method for sparse precision matrix  #
#             estimation                                                           #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 23th, 2013                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#


camel.tiger.clime.mfista <- function(Sigma, d, maxdf, mu, lambda, shrink, prec, max.ite){
  
  d_sq = d^2
  Y = diag(d)
#   lambda = lambda-shrink*prec
  nlambda = length(lambda)
  L = eigen(Sigma)$values[1]^2
  icov = array(0,dim=c(d,d,nlambda))
  ite.ext = rep(0,d*nlambda)
  obj = array(0,dim=c(max.ite,nlambda))
  runt = array(0,dim=c(max.ite,nlambda))
  x = array(0,dim=c(d,maxdf,nlambda))
  col_cnz = rep(0,d+1)
  row_idx = rep(0,d*maxdf*nlambda)
  begt=Sys.time()
  str=.C("tiger_clime_mfista", as.double(Y), as.double(Sigma), as.double(icov), 
         as.integer(d), as.double(mu), as.integer(ite.ext), as.double(lambda), 
         as.integer(nlambda), as.integer(max.ite), as.double(prec), as.double(L), 
         as.double(x), as.integer(col_cnz), as.integer(row_idx),PACKAGE="camel")
  runt1=Sys.time()-begt
  ite.ext = matrix(unlist(str[6]), byrow = FALSE, ncol = nlambda)
  obj = 0
  icov_list = vector("list", nlambda)
  icov_list1 = vector("list", nlambda)
  for(i in 1:nlambda){
    icov_i = matrix(unlist(str[3])[((i-1)*d_sq+1):(i*d_sq)], byrow = FALSE, ncol = d)
    icov_list1[[i]] = icov_i
    icov_list[[i]] = icov_i*(abs(icov_i)<=abs(t(icov_i)))+t(icov_i)*(abs(t(icov_i))<abs(icov_i))
    obj[i] = sum(abs(icov_i))
  }
  x = unlist(str[12])
  col_cnz = unlist(str[13])
  row_idx = unlist(str[14])
  return(list(icov=icov_list, icov1=icov_list1,ite=ite.ext, obj=obj,runt=runt1,
              x=x, col_cnz=col_cnz, row_idx=row_idx))
}
