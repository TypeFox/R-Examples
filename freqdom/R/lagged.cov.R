#' For a given multivariate stationary time series estimates a covarianve matrix
#' \eqn{C_{XY}^k = Cov(X_k,Y_0)} using the formula
#' \deqn{\hat C_{XY}^k =  \frac{1}{n} \sum_{i=1}^{n-k} X_{k+i} Y_k'. }
#'
#' @title Compute cross covariance with a given lag
#' @param X first process 
#' @param Y second process, if null then autocovariance of X is computed
#' @param lag the lag that we are interested in
#' @return Covariance matrix 
#' @export
#' @examples
#' X = rar(100)
#' Y = rar(100)
#' lagged.cov(X,Y)
lagged.cov = function(X,Y=NULL,lag=0){
  if (base::is.null(Y))
		Y = X

  if (base::dim(X)[1] != base::dim(Y)[1])
    stop("Number of observations must be equal")
  if (!base::is.matrix(X) || !base::is.matrix(Y))
    stop("X and Y must be matrices")
  
	n = base::dim(X)[1]
	h = base::abs(lag)
	
  if (n - 1 <= h)
	  base::stop(base::paste("Too little observations to compute lagged covariance with lag",h))
	
  M = base::t(X[1:(n-h),]) %*% (Y[1:(n-h)+h,])/(n)
	if (lag < 0){
	  M = base::t(Y[1:(n-h),]) %*% (X[1:(n-h)+h,])/(n)
	  M = base::t(M)
	}
	M
}

