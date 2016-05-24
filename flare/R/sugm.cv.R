#----------------------------------------------------------------------------------#
# Package: flare                                                                   #
# sugm.cv(): Cross validation for regularization parameter                         #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Dec 2nd 2013                                                               #
# Version: 1.1.0                                                                   #
#----------------------------------------------------------------------------------#

sugm.cv <- function(obj, loss=c("likelihood", "tracel2"), fold=5) {
  x = obj$data
  if (is.null(x)) stop("No data matrix in sugm object.  Use data matrix instead of sample covariance matrix for computing sugm!")
  n = nrow(x)
  p = ncol(x)
  
  part.list = part.cv(n, fold)

  
  lossname = match.arg(loss, c("likelihood", "tracel2"))
  lossname = paste("sugm", lossname, sep=".")
  lossfun = match.fun(lossname)
  
  loss.re = matrix(0, nrow = fold, ncol = obj$nlambda)
  scalar = 1-1/nrow(x[part.list$testMat[,1],])
  for (i in 1:fold) {
    x.train = x[part.list$trainMat[,i],]
    sugm.cv = sugm(x.train, lambda=obj$lambda, method=obj$method,sym=obj$sym,verbose=obj$verbose,
                     standardize=obj$standardize)
    x.test = x[part.list$testMat[,i],]
    ntest = nrow(x.test)
    for (j in 1:obj$nlambda) {
      loss.re[i,j] = loss.re[i,j]  + lossfun(cov(x.test)*scalar,  sugm.cv$icov[[j]])
    }
  }
  
  loss.mean = apply(loss.re, 2, mean)
  loss.sd = apply(loss.re, 2, sd)
  
  opt.idx = which.min(loss.mean)
  lambda.opt = obj$lambda[opt.idx]
  
  outlist = list(lambda.opt=lambda.opt, opt.idx=opt.idx,loss=lossname, loss.mean=loss.mean, loss.sd = loss.sd)
  class(outlist) = c("sugm.cv")
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

sugm.likelihood <- function(Sigma, Omega) {
  ot = as.numeric(unlist(determinant(Omega)))
  if (ot[2]<=0) warning("Precision matrix estimate is not positive definite!")
  tmp = (sum(diag(Sigma%*%Omega))  - ot[1])
  if(is.finite(tmp)) {
    return(tmp)
  } else {
    return(Inf)
  }
}

sugm.tracel2 <- function(Sigma, Omega) {
  return(sum(diag(  (Sigma%*%Omega -  diag(1,  dim(Omega)[1] ) )^2 )) )
}
