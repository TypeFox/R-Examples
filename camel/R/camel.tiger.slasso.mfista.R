#----------------------------------------------------------------------------------#
# Package: camel                                                                   #
# camel.tiger.slasso.mfista(): Scaled Lasso method for sparse precision matrix      #
#             estimation                                                           #
# Authors: Xingguo Li                                                              #
# Emails: <xingguo.leo@gmail.com>                                                  #
# Date: Aug 23th, 2013                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#


camel.tiger.slasso.mfista <- function(data, n, d, maxdf, mu, lambda, shrink, prec, max.ite){
  
  nlambda = length(lambda)
  lambda = lambda-shrink*prec
  X1=data
  d_sq=d^2
  X1 = X1 - matrix(rep(colMeans(X1),n), nrow=n, byrow=TRUE)
  Gamma=diag(diag(t(X1)%*%X1/n))
  Q = diag(1/sqrt(diag(Gamma)))
  Omega=array(0,dim=c(d,d,nlambda))
  Z=X1%*%(diag(1/sqrt(diag(Gamma))))
  ZZ = t(Z)%*%Z
  L=eigen(ZZ)$values[1]
  
  icov = array(0,dim=c(d,d,nlambda))
  ite_ext = rep(0,d*nlambda)
  x = array(0,dim=c(d,maxdf,nlambda))
  col_cnz = rep(0,d+1)
  row_idx = rep(0,d*maxdf*nlambda)
  icov_list = vector("list", nlambda)
  icov_list1 = vector("list", nlambda)
  obj = array(0,dim=c(max.ite,nlambda))
  runt = array(0,dim=c(max.ite,nlambda))
  str=.C("tiger_slasso_mfista", as.double(Z), as.double(icov), as.double(x), as.integer(d), as.integer(n), 
         as.double(mu), as.integer(ite_ext), as.double(lambda), as.integer(nlambda), as.integer(max.ite), 
         as.integer(col_cnz),as.integer(row_idx), as.double(prec), as.double(L), 
         as.double(obj),as.double(runt),PACKAGE="camel")
  for(i in 1:nlambda){
    icov_i = Q%*%(matrix(unlist(str[2])[((i-1)*d_sq+1):(i*d_sq)], byrow = FALSE, ncol = d))%*%Q
    icov_list1[[i]] = icov_i
    icov_list[[i]] = icov_i*(abs(icov_i)<=abs(t(icov_i)))+t(icov_i)*(abs(t(icov_i))<abs(icov_i))
  }
  ite_ext = matrix(unlist(str[7]), byrow = FALSE, ncol = nlambda)
  x = unlist(str[3])
  col_cnz = unlist(str[11])
  row_idx = unlist(str[12])
  return(list(icov=icov_list, icov1=icov_list1, ite=ite_ext, x=x, col_cnz=col_cnz, row_idx=row_idx))
}
