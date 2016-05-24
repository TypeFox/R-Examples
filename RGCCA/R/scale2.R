#' Standardization (to zero means and unit variances) of matrix-like objects.
#' @param A A numeric matrix.
#' @param center A logical value. If center = TRUE, each column is translated to have zero mean (default: TRUE).
#' @param scale A logical value. If scale = TRUE, each column is transformed to have unit variance (default = TRUE).
#' @param bias Logical value for biaised (\eqn{1/n}) or unbiaised (\eqn{1/(n-1)}) estimator of the var/cov (default = TRUE).
#' @return \item{A}{The centered and/or scaled matrix. The centering and scaling values (if any) are returned as attributes "scaled:center" and "scaled:scale".}
#' @title Scaling and Centering of Matrix-like Objects
#' @export scale2

scale2<-function (A, center = TRUE, scale = TRUE, bias = TRUE) 
{
  if (center == TRUE & scale == TRUE) {
    A = scale(A, center = TRUE, scale = FALSE)
    std = sqrt(apply(A, 2, function(x) cov2(x, bias = bias)))
    if (any(std==0)) {
      sprintf("there were %d constant variables",sum(std==0))
      std[std==0]=1
    }
    A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    return(A)
  }
  if (center == TRUE & scale == FALSE) {
    A = scale(A, center = TRUE, scale = FALSE)
    return(A)
  }
  if (center == FALSE & scale == TRUE) {
    std = apply(A, 2, function(x) cov2(x, bias = bias))
    A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    return(A)
  }
}
