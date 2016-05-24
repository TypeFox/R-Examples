
#########
#### Extracting correlations from a covariance matrix
#########
#' Extracting correlations from a covariance matrix
#'
#' @param mat A covariance matrix.
#' @return The correlation matrix embedded in \code{mat}.
#' @examples
#' # 2 dimensional case
#' d <- 2
#' tmp <- matrix(rnorm(d^2), d, d)
#' mcov <- tcrossprod(tmp, tmp)
#' 
#' # Covariance matrix
#' mcov
#' 
#' # Correlation matrix
#' extractCorr(mcov)
#' @export

extractCorr <- function(mat)
{
  diag(sqrt( diag(mat) )^-1, nrow(mat)) %*% mat %*% diag(sqrt( diag(mat) )^-1, nrow(mat))  
}
