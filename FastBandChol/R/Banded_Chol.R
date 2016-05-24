#' Computes estimate of covariance matrix by banding the Cholesky factor

#' @param X an n by p data matrix
#' @param bandwidth an integer less than n-1
#' @param centered TRUE/FALSE; is data matrix centered?

#' @return est estimated covariance matrix



banded.chol = function(X, bandwidth, centered=FALSE){
  n = nrow(X)
  p = ncol(X)
  

  if(bandwidth > n-1){
    stop("bandwidth larger than n-1; solution undefined")
  }
  
  if(bandwidth > p-1){
    stop("bandwidth larger than p-1; consider smaller bandwidth")
  }
  

  if(!centered){
    mx = apply(X, 2, mean)
      x.c = X - tcrossprod(rep(1, n), mx)
    } else{
      x.c = X
    }

    R = BandCholcpp(x.c, bandwidth)$R
    est = n^(-1)*crossprod(R)

    return(list("est" = est))
}
