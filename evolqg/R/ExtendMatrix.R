#' Control Inverse matrix noise with Extension
#'
#' Calculates the extented covariance matrix estimation as described in Marroig et al. 2012
#'
#' @param cov.matrix Covariance matrix
#' @param var.cut.off Cut off for second derivative variance. Ignored if ret.dim is passed.
#' @param ret.dim Number of retained eigen values
#' @return Extended covariance matrix and second derivative variance
#' @references Marroig, G., Melo, D. A. R., and Garcia, G. (2012). Modularity, noise, and natural selection. Evolution; international journal of organic evolution, 66(5), 1506-24. doi:10.1111/j.1558-5646.2011.01555.x
#' @author Diogo Melo
#' @note Covariance matrix being extended should be larger then 10x10
#' @export
#' @examples
#' cov.matrix = RandomMatrix(11, 1, 1, 100)
#' ext.matrix = ExtendMatrix(cov.matrix, var.cut.off = 1e-6)
#' ext.matrix = ExtendMatrix(cov.matrix, ret.dim = 6)
#' @keywords extension
#' @keywords covariancematrix
ExtendMatrix <- function(cov.matrix, var.cut.off = 1e-4, ret.dim = NULL) {
  p = dim(cov.matrix)[1]
  if(p<10)
    warning("matrix is too small")
  eigen.cov.matrix = eigen(cov.matrix)
  eVal = eigen.cov.matrix$values
  eVec = eigen.cov.matrix$vectors
  
  grad = array(dim=c(p-2))
  tr.cov.matrix = sum(eVal)
  for (i in 1:(p-2))
    grad[i] = abs(eVal[i]/tr.cov.matrix - 2*(eVal[i+1]/tr.cov.matrix) + eVal[i+2]/tr.cov.matrix)
  var.grad = array(dim=c(p-6))
  for (i in 1:(p-6)){
    var.grad[i] = var(grad[i:(i+4)])
  }
  if(is.null(ret.dim)){
    ret.dim = which(var.grad < var.cut.off)[1]
  }
  eVal[eVal < eVal[ret.dim]] = eVal[ret.dim]
  extended.cov.matrix = eVec%*%diag(eVal)%*%t(eVec)
  colnames(extended.cov.matrix) = colnames(cov.matrix)
  rownames(extended.cov.matrix) = rownames(cov.matrix)
  return(list(ExtMat = extended.cov.matrix, var.grad = var.grad, eVals = eVal))
}
