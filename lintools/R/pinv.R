
#' Moore-Penrose pseudoinverse
#'
#' Compute the pseudoinverse of a matrix using the
#' SVD-construction 
#' 
#' 
#' @section Details:
#' 
#' The Moore-Penrose pseudoinverse (sometimes called the generalized inverse) \eqn{\boldsymbol{A}^+} of a matrix \eqn{\boldsymbol{A}}
#' has the property that \eqn{\boldsymbol{A}^+\boldsymbol{AA}^+ = \boldsymbol{A}}. It can be constructed as follows.
#' 
#' \itemize{
#' \item{Compute the singular value decomposition \eqn{\boldsymbol{A} = \boldsymbol{UDV}^T}}
#' \item{Replace diagonal elements in \eqn{\boldsymbol{D}} of which the absolute values are larger than some limit \code{eps} with their reciprocal values}
#' \item{Compute \eqn{\boldsymbol{A}^+ = \boldsymbol{UDV}^T}}
#' } 
#' 
#'
#' @param A [numeric] matrix
#' @param eps [numeric] tolerance for determining zero singular values
#'
#' @references
#' 
#' S Lipshutz and M Lipson (2009) Linear Algebra. In: Schuam's outlines. McGraw-Hill
#' 
#' @examples 
# example from appendix of Lipschutz and Lipson (2009)
#' A <- matrix(c(
#'  1,  1, -1,  2,
#'  2,  2, -1,  3,
#'  -1, -1,  2, -3
#' ),byrow=TRUE,nrow=3)
#' # multiply by 55 to retrieve whole numbers
#' pinv(A) * 55
#'
#' @export
pinv <- function(A, eps=1e-8){
  L <- svd(A)
  d <- L$d
  i <- abs(d) > eps
  d[i] <- 1/d[i]
  L$v %*% diag(d, nrow=length(d)) %*% t(L$u)
}















  
