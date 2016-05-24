#' Calculate the conditional mean & variance of a random vector following a multivariate normal distribution
#'
#' Calculate the conditional mean & variance of a sub-vector
#' of a random vector given the rest elements, where the random vector follows a multivariate normal distribution.
#' This function calculates the conditional mean & variance of a multivariate normal distribution,
#' with (\code{meanVec}, \code{varMat}), conditional on y at the indices \code{conditionalIndex}.
#'
#' @param y the vector of the known values to be conditioned on.
#' @param conditionalIndex the index to be conditioned on.
#' @param meanVec the mean vector the joint multivariate normal distribution.
#' @param varMat the variance-covariance matrix of the joint multivariate normal distribution.
#' @return a list consisting of 'mean' and 'var' representing the conditional mean and variance respectively.
#' @keywords internal
# @examples
# conditionalDistMvnorm(c(-0.5,0.5), c(2,4), c(1,2,3,4),matrix(
# c(1,0.3,0.2,0.1, 0.3,1,-0.1,0.3,0.2,-0.1,1,0.1,0.1,0.3,0.1,1),
# nrow=4,ncol=4,byrow=TRUE))

conditionalDistMvnorm <- function(y, conditionalIndex, meanVec, varMat)
{
	sigma11 <- varMat[-conditionalIndex,-conditionalIndex]
	sigma22 <- varMat[ conditionalIndex, conditionalIndex]
	sigma12 <- varMat[-conditionalIndex, conditionalIndex]
	sigma21 <- varMat[ conditionalIndex,-conditionalIndex]
	invSigma22 <- solve(sigma22)
	mNew <- meanVec[-conditionalIndex] + as.vector(sigma12 %*% invSigma22 %*%(y - meanVec[conditionalIndex]))
	vNew <- sigma11 - sigma12 %*% invSigma22 %*% sigma21

	return (list('mean'=mNew,'var' = vNew))
}

isStationaryAR <- function(arPrmtr)
{
  roots <- polyroot(c(1,-arPrmtr))
	#print(roots)
	if(any(abs(roots)<=1))
    FALSE
  else
    TRUE
}



#' Compute the covariance matrix of some observations of the AR model
#'
#' Compute the covariance matrix of \eqn{(\eta_t,...,\eta_{t-lag})} for an AR model.
#' @param arPrmtr the parameter of the AR model (with no intercept term).
#' @param sigma the innovation standard deviation.
#' @param lag the number of lags to be computed, including lag zero.
#' @return the covariance matrix.
#' @keywords internal
computeCovAR <- function(arPrmtr, sigma, lag=length(arPrmtr)+1)
{
	roots <- polyroot(c(1,-arPrmtr))
	#print(roots)
	if(any(abs(roots)<=1))
		warning(" arPrmtr is not stationary",paste(arPrmtr,sep=',')," Roots:", roots)
	val <- stats::ARMAacf(ar=arPrmtr, lag.max=max(length(arPrmtr),lag))
	#print(val)
	val <- as.vector(val)

	mat <- matrix(nrow=lag,ncol=lag)
	for(i in 1:lag)
	{
		mat[i,i] <- 1
		if(i > 1)
    {
			for(j in 1:(i-1))
			{
				mat[i,j] <- val[abs(i-j)+1]
				mat[j,i] <- mat[i,j]
			}
		}
	}

	v <- sigma^2/(1-arPrmtr%*%val[2:(length(arPrmtr)+1)])
	mat <- v[1,1]*mat
	return(mat)
}
