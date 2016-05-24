#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# sugm.clime.ladm.scr(): Coordinate descent method for sparse precision matrix     #
#                        estimation                                                #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Jul 8th, 2014                                                              #
# Version: 1.4.0                                                                   #
#----------------------------------------------------------------------------------#

sugm.clime.ladm.scr <- function(Sigma, lambda, nlambda, n, d, maxdf, rho, shrink, prec, max.ite, verbose){
  
  if(verbose==TRUE)
    cat("The Constrained L1 Minimization for Sparse Precision Matrix Estimation.\n")
  d.sq = d^2
  lambda = lambda-shrink*prec
  idx.scr = apply(Sigma,2,order,decreasing=TRUE)
  num.scr = d
  if(d>=n){
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
  ite.int = matrix(0,nrow=d,ncol=nlambda)
  ite.int1 = matrix(0,nrow=d,ncol=nlambda)
  ite.int2 = matrix(0,nrow=d,ncol=nlambda)
  x = rep(0,d*maxdf*nlambda)
  col.cnz = rep(0,d+1)
  row.idx = rep(0,d*maxdf*nlambda)
  icov.list1 = vector("list", nlambda)
  for(i in 1:nlambda){
    icov.list1[[i]] = matrix(0,d,d)
  }
  for(j in 1:d){
    idx.scr0 = idx.scr[,j]
    idx.scr1 = idx.scr0[1:num.scr1]
    idx.scr2 = idx.scr0[1:num.scr2]
    S.order = Sigma[idx.scr0,idx.scr0]
    SS.order = crossprod(S.order)
    gamma = max(colSums(abs(SS.order)))
    icov0 = rep(0,d*nlambda)
    ite0.int = rep(0,nlambda)
    ite0.int1 = rep(0,nlambda)
    ite0.int2 = rep(0,nlambda)
    x0 = rep(0,maxdf*nlambda)
    col.cnz0 = 0
    row.idx0 = rep(0,maxdf*nlambda)
    str=.C("sugm_clime_ladm_scr", as.double(Sigma), as.double(SS.order), as.double(icov0), as.double(x0), 
           as.integer(d), as.double(gamma), as.double(lambda), as.integer(nlambda), 
           as.double(rho), as.integer(col.cnz0), as.integer(row.idx0), 
           as.integer(ite0.int), as.integer(ite0.int1), as.integer(ite0.int2), 
           as.integer(num.scr1), as.integer(num.scr2), 
           as.integer(idx.scr0), as.integer(idx.scr1), as.integer(idx.scr2), 
           as.integer(max.ite), as.double(prec), as.integer(j), PACKAGE="flare")
    icov = matrix(unlist(str[3]), byrow = FALSE, ncol = nlambda)
    for(i in 1:nlambda){
      icov.list1[[i]][,j] = icov[,i]
    }
    cnt = unlist(str[10])
    col.cnz[j+1] = cnt+col.cnz[j]
    x[(col.cnz[j]+1):col.cnz[j+1]] = unlist(str[4])[1:cnt]
    row.idx[(col.cnz[j]+1):col.cnz[j+1]] = unlist(str[11])[1:cnt]
    ite.int[j,] = unlist(str[12])
    ite.int1[j,] = unlist(str[13])
    ite.int2[j,] = unlist(str[14])
  }
  icov.list = vector("list", nlambda)
  for(i in 1:nlambda){
    icov.i = icov.list1[[i]]
    icov.list[[i]] = icov.i*(abs(icov.i)<=abs(t(icov.i)))+t(icov.i)*(abs(t(icov.i))<abs(icov.i))
  }
  ite = list()
  ite[[1]] = ite.int1
  ite[[2]] = ite.int2
  ite[[3]] = ite.int
  return(list(icov=icov.list, icov1=icov.list1,ite=ite, x=x, col.cnz=col.cnz, row.idx=row.idx))
}
