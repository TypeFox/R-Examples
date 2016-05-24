#' Selects tuning parameter for banded sample covariance matrix by k fold cross-validation

#' @param X an n by p data matrix
#' @param bandwidth a vector of candidate bandwidths
#' @param folds number of folds for CV; default is 3
#' @param est.eval TRUE/FALSE should the solution at the minizing tuning paramter be computed
#' @param Frob TRUE/FALSE; which type of norm used for evaluting CV error; FALSE sets to operator norm

#' @return bandwidth.min the bandwidth minimizing cv error
#' @return est the banded estimate


banded.sample.cv = function(X, bandwidth, folds=3, est.eval=TRUE, Frob=TRUE){
  if(folds!=3){
    f = folds
  } else {
    f = 3
  }

  n = nrow(X); p = ncol(X);
  nk = length(bandwidth)

  fold = sample(rep(1:f, length=n))
  cv.index = split(1:n, fold)

  cv.func = function(cv.index){
    m = length(cv.index)
    x.test = X[cv.index,]
    mx = apply(x.test, 2, mean)
    x.c.test = x.test - tcrossprod(rep(1, m), mx)
    S.test = m^(-1)*crossprod(x.c.test)
    x.train = X[-cv.index,]
    err = numeric(length=nk)

    for(i in 1:nk){
      est = banded.sample(x.train, bandwidth=bandwidth[i], centered=FALSE)$est
      if(!Frob){
        err[i] = abs(max(eigen(S.test - est)$val))
      } else {
        err[i] = norm(S.test - est, type="F")^2
      }
    }

    return(err)
  }

  result = lapply(cv.index, cv.func)
  result = colSums(do.call(rbind, result))
  bandwidth.ind = which.min(result)
  bandwidth.min = bandwidth[bandwidth.ind]

  if(!est.eval){
    return(list("bandwidth.min" = bandwidth.min))
  } else {
    est = banded.sample(X, bandwidth=bandwidth.min, centered=FALSE)$est
    return(list("est" = est, "bandwidth.min" = bandwidth.min))
  }
}



