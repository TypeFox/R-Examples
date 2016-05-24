#' Compute
#' \deqn{Y_t = \sum_{k=-q}^p A_k X_{t-k}}
#' where \eqn{X_t} is a stationary multivariate time series and \eqn{(A_k)_{-q \leq k \leq p}} is a filter.
#'
#' @title Generate a linear process from 
#' @param X process process
#' @param A time-domain operator series
#' @param noise function taking dimension D and returning D-dimentional vector
#' @return Multivariate linear process
#' @importFrom stats rnorm
#' @export
#' @seealso \code{\link{speclagreg}}
#' @examples
#' d = 2
#' n = 100
#' X = rar(n,d=d)
#' 
#' OP = array(0,c(d,d,2))
#' OP[,,1] = 2 * diag(d:1)/d
#' OP[,,2] = 1.5 * diag(d:1)/d
#' A = timedom(OP, 0:1)
#' 
#' Y = linproc(X,A,noise=rnorm)
linproc = function(X, A, noise=NULL){
  if (!is.matrix(X))
    stop("X must be matrix")
  if (!is.timedom(A))
    stop("A must be a time domain operator")
  
  # if no noise then gaussian noise
	if (is.null(noise))
		noise = function(n){ rep(0,n) }
  WN = c()
  for (i in 1:dim(X)[1])
    WN = rbind(WN,noise(dim(A$operators)[2]))
  
  # concolution + noise
  X = A %c% X
  X + WN
}
