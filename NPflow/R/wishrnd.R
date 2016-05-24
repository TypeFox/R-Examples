#' Sample from a Wishart distribution
#'
#' For internal use only.
#'
#'@param n degrees of freedom
#'
#'@param Sigma scale parameter (matrix)
#'
#'@keywords internal
#'
#'@importFrom stats rnorm
#'
#'@export wishrnd

wishrnd <- function(n, Sigma){

  p <- dim(Sigma)[1]
  p2 <- dim(Sigma)[2]

  if(p!=p2){
    stop('Error : Matrix not square\n')
  }

  x=matrix(rnorm(n=n*p), nrow=n, ncol=p)%*%chol(Sigma)
  W=t(x)%*%x

  return(W)
}
