#----------------------------------------------------------------------------------#
# Package: camel                                                                   #
# camel.tiger.cv(): Cross validation for regularization parameter                  #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Aug 23th, 2013                                                             #
# Version: 0.1.0                                                                   #
#----------------------------------------------------------------------------------#

camel.tiger.cv <- function(obj, loss=c("likelihood", "tracel2"), fold=5) {
  x = obj$data
  if (is.null(x)) stop("No data matrix in tiger object.  Use data matrix instead of sample covariance matrix for computing tiger!")
  n = nrow(x)
  p = ncol(x)
  
  part_list = part.cv(n, fold)

  
  lossname = match.arg(loss, c("likelihood", "tracel2"))
  lossname = paste("tiger", lossname, sep=".")
  lossfun = match.fun(lossname)
  
  loss_re = matrix(0, nrow = fold, ncol = obj$nlambda)
  scalar = 1-1/nrow(x[part_list$testMat[,1],])
  for (i in 1:fold) {
    x_train = x[part_list$trainMat[,i],]
    tiger_cv = camel.tiger(x_train, lambda=obj$lambda, method=obj$method,sym=obj$sym,verbose=obj$verbose,
                     standardize=obj$standardize,correlation=obj$correlation)
    x_test = x[part_list$testMat[,i],]
    ntest = nrow(x_test)
    for (j in 1:obj$nlambda) {
      loss_re[i,j] = loss_re[i,j]  + lossfun(cov(x_test)*scalar,  tiger_cv$icov[[j]])
    }
  }
  
  loss_mean = apply(loss_re, 2, mean)
  loss_sd = apply(loss_re, 2, sd)
  
  opt_idx = which.min(loss_mean)
  lambda_opt = obj$lambda[opt_idx]
  
  outlist = list(lambda_opt=lambda_opt, opt_idx=opt_idx,loss=lossname, loss_mean=loss_mean, loss_sd = loss_sd)
  class(outlist) = c("tiger.cv")
  return(outlist)
}

part.cv <- function(n, fold) {
  
  ntest = floor(n/fold)
  ntrain = n-ntest
  
  ind = sample(n)
  
  trainMat = matrix(NA, nrow=ntrain, ncol=fold)
  testMat = matrix(NA, nrow=ntest, ncol=fold)
  
  nn = 1:n
  
  for (j in 1:fold) {
    sel = ((j-1)*ntest+1):(j*ntest)
    testMat[,j] = ind[sel ]
    sel2 =nn[ !(nn %in% sel) ]
    trainMat[,j] = ind[sel2]
  }
  
  return(list(trainMat=trainMat, testMat=testMat))
}

tiger.likelihood <- function(Sigma, Omega) {
  ot = as.numeric(unlist(determinant(Omega)))
  if (ot[2]<=0) warning("Precision matrix estimate is not positive definite!")
  tmp = (sum(diag(Sigma%*%Omega))  - ot[1])
  if(is.finite(tmp)) {
    return(tmp)
  } else {
    return(Inf)
  }
}

tiger.tracel2 <- function(Sigma, Omega) {
  return(sum(diag(  (Sigma%*%Omega -  diag(1,  dim(Omega)[1] ) )^2 )) )
}
