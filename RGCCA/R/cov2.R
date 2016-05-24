#' cov2() is similar to cov() but has an additional argument. The denominator \eqn{n} (bias = TRUE) 
#' can be used (instead of \eqn{n-1}) to give a biased estimator of the (co)variance.
#' @param x A numeric vector, matrix or data.frame.
#' @param y A numeric vector, matrix or data.frame.
#' @param bias A logical value. If bias = TRUE, \eqn{n} is used to give a biased estimator of the (co)variance. 
#' If bias = FALSE, \eqn{n-1} is used (default: TRUE).
#' @return \item{C}{Estimation of the variance (resp. covariance) of x (resp. x and y).}
#' @title Variance and Covariance (Matrices)
#' @export cov2

cov2 = function (x, y = NULL, bias = TRUE) 
{
  n = NROW(x)
  
  if (is.null(y)) {
    x = as.matrix(x)
    if(bias){
      C = ((n - 1)/n) * cov(x, use = "pairwise.complete.obs")  
    }
    else{
      C = cov(x, use = "pairwise.complete.obs") 
    }
  }
  
  else{
    if(bias){
      C = ((n - 1)/n) * cov(x, y, use = "pairwise.complete.obs")  
    }
    else{
      C = cov(x, y, use = "pairwise.complete.obs") 
    }
  }
return(C)
}
