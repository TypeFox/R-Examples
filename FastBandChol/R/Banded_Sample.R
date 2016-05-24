#' Estimates a covariance matrix by banding the sample covariance matrix

#' @param X an n times p data matrix
#' @param bandwidth nonnegative integer; less than n-1
#' @param centered TRUE/FALSE is data matrix centered?

#' @return est the banded estimate


banded.sample = function(X, bandwidth, centered=FALSE){
  n = nrow(X)
  p = ncol(X)

  if(!centered){
    mx = apply(X, 2, mean)
    x.c = X - tcrossprod(rep(1, n), mx)
  } else{
    x.c = X
  }

  S = n^(-1)*crossprod(x.c)

  for(l in 1:p){
    for(r in 1:p){
      if(abs(l-r) > bandwidth){
          S[l,r] = 0
        } else{
          }
      }
    }

  return(list("est" = S))
}
