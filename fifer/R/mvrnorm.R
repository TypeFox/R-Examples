##' Randomly Generate Multivariate Normal Data
##'
##' This function will randomly generate correlated multivariate normal data with specified means
##' and covariances (or correlations). The user also has the flexibility to generate data with 
##' a randomly selected correlation matrix using the \code{\link{random.correlation}} function.
##'	
##' \code{mv.norm} generates correlated multivariate normal data using a choleski decomposition. If the
##' user does not specify \code{Sigma}, a random correlation matrix will be generated. Also, if means are
##' not specified, the function will default to means of zero. 
##' @param n The sample size of the randomly generated dataset
##' @param vars An integer indicating the number of variables. Ignored unless no \code{Sigma} is supplied.
##' @param mu A vector of means that has the same length as the number of rows/columns of Sigma. Defaults to 
##' a vector of zeroes. 
##' @param Sigma A positive definate matrix. If \code{NULL}, the user must specify \code{vars}. 
##' @param names Optional. A vector of strings indicating the variable names. 
##' @aliases mv.Rnorm mvrnorm 
##' @seealso \code{\link{random.correlation}} \code{\link{cor2cov}}
##' @return a nxp matrix of pseuodo-random values. 
##' @author Dustin Fife
##' @export
##' @examples
##' set.seed(2)
##' ## generate data with correlation of .6
##' d = mv.rnorm(n=1000, Sigma=matrix(c(1, .6, .6, 1), 2), names=c("x", "y"))
##' head(d); cor(d)
##' ## generate data with a random correlation
##' d = mv.rnorm(n=1000, vars=4, names=letters[1:4])
##' head(d); cor(d)
##' ## generate non-scaled data
##' ms = c(100, 10, 5, 0) ### specify means
##' Sigma = matrix(c(1, .6, .5, .4,
##' 			.6, 1, .3, .2,
##' 			.5, .3, 1, .1,
##' 			.4, .2, .1, 1), 4)
##' ## convert Sigma to covariance matrix
##' Sigma = cor2cov(Sigma, sd=c(15, 3, 2, 1))
##' ## generate the data
##' d = mv.rnorm(n=1000, mu=ms, Sigma=Sigma, names=letters[1:4])
##' head(d); cor(d)
mv.rnorm = function(n=1, vars=NULL, mu=NULL, Sigma=NULL, names=NULL){

	#### if they supply vars and no Sigma, randomly generate correlation matrix
	if (is.null(Sigma)){
		if (is.null(vars)){
			stop("If Sigma is set to NULL, vars must not be NULL.")
		} else {
			Sigma = random.correlation(vars)	
		}
	} 
	
	p = ncol(Sigma)
	
	if (is.null(mu)){
		mu = rep(0, times=p)
	}
	
	#### bark if the dimensions do not match
	if (dim(Sigma)[1] != length(mu) | ncol(Sigma) != nrow(Sigma)){
		stop("Either Sigma is not symmetric or mu and Sigma differ in dimensions.")
	}

	#### generate uncorrelated data
	X = matrix(rnorm(p*n), n)
	
	#### perform choleski decomp
	a = chol(Sigma)
	X = X%*%a
	
	#### scale it to have the right means
	sd = diag(Sigma)
	X = mu + sd*scale(X)
	
	#### name the columns (if user specifies them)
	if (!is.null(names)){
		X = data.frame(X)
		names(X) = names
	}	

	
	#### return the matrix
	return(X)
}