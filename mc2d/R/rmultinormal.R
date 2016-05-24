#<<BEGIN>>
rmultinormal <- function(n, mean, sigma , method=c("eigen", "svd", "chol"))
#TITLE The Vectorized Multivariate Random Deviates
#NAME multinormal
#KEYWORDS distribution
#DESCRIPTION
#This function is the vectorized version of the \samp{rmvnorm} from the \samp{mvtnorm} library.
#It provides a random number generator for the multivariate normal distribution 
#with varying vectors of means and varying covariance matrixes.
#INPUTS
#{x}<<Vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.>> 
#{n}<<Number of observations. If \samp{length(n) > 1}, the length is taken to be the number required.>>
#{mean}<<Vector or matrix of means. If a matrix, each row is taken to be a quantile.
#Default is a vector of 0 of convenient length.>>
#{sigma}<<Covariance vector corresponding to the coercion of the covariance matrix into a vector (if unique for all \samp{n} or \samp{x}) 
#or array of covariance vectors (if varying according to \samp{n} or \samp{x}).
#default is a diagonal matrix of convenient sizee.>>
#{method}<<Matrix decomposition used to determine the matrix root of sigma, possible methods are
#eigenvalue decomposition ("eigen", default), singular value decomposition ("svd"), and Cholesky decomposition ("chol").>>
#{log}<<Logical; if \samp{TRUE}, densities d are given as log(d).>> 
#DETAILS
#\samp{rmvnorm(n, m, s)} is equivalent to \samp{rmultinormal(n, m, as.vector(s))}.
#\samp{dmvnorm(x, m, s)} is equivalent to \samp{dmultinormal(x, m, as.vector(s))}.
#
#If \samp{mean} and/or \samp{sigma} is a matrix, 
#the first random deviate will use the first row of \samp{mean} and/or \samp{sigma}, the second random
#deviate will use the second row of \samp{mean} and/or \samp{sigma}, ...
#recycling being permitted by raw.
#If \samp{mean} is a vector of length \samp{l} or is a matrix with \samp{l} columns, \samp{sigma}
#should be a vector of length \samp{l x l} or a matrix of number of \samp{l x 2} columns. 
#NOTE
#The use of a varying sigma may be very time consumming.
#EXAMPLE
#
### including equivalence with dmvnorm
### mean and sigma as vectors
#(mean <- c(10,0))
#(sigma <- matrix(c(1,2,2,10), ncol=2))
#sigma <- as.vector(sigma)
#(x <- matrix(c(9,8,1,-1), ncol=2))
#round(rmultinormal(10,mean,sigma))
#dmultinormal(x,mean,sigma)               
### Eq
#dmvnorm(x,mean,matrix(sigma, ncol=2))               
#
### mean as matrix
#(mean <- matrix(c(10,0,0,10),ncol=2))
#round(rmultinormal(10,mean,sigma))
#dmultinormal(x,mean,sigma)
### Eq
#dmvnorm(x[1,],mean[1,],matrix(sigma, ncol=2))               
#dmvnorm(x[2,],mean[2,],matrix(sigma, ncol=2))               
#
### sigma as matrix
#(mean <- c(10,0))
#(sigma <- matrix(c(1,2,2,10,10,2,2,1), nrow=2, byrow=TRUE))
#round(rmultinormal(10,mean,sigma))
#dmultinormal(x,mean,sigma)               
### Eq
#dmvnorm(x[1,],mean,matrix(sigma[1,],ncol=2))               
#dmvnorm(x[2,],mean,matrix(sigma[2,],ncol=2))               
#
### mean and sigma as matrix
#(mean <- matrix(c(10,0,0,10),ncol=2))
#(sigma <- matrix(c(1,2,2,10,10,2,2,1), nrow=2, byrow=TRUE))
#round(rmultinormal(10,mean,sigma))
#dmultinormal(x,mean,sigma)               
### Eq
#dmvnorm(x[1,],mean[1,],matrix(sigma[1,],ncol=2))               
#dmvnorm(x[2,],mean[2,],matrix(sigma[2,],ncol=2))               
#
#(mean <- c(10,0))
#(sigma <- matrix(c(1,2,2,10,10,2,2,1), nrow=2, byrow=TRUE))
#x <- rmultinormal(1000,mean,sigma)
#plot(x)
#--------------------------------------------
{
  if(length(n) == 0) return(n)
  if(length(n) > 1) n <- length(n)
  if (missing(mean)){
	if(is.vector(sigma))  mean <- rep(0, length = sqrt(length(sigma)))
	else  mean <- rep(0, length = sqrt(ncol(sigma)))}
  if (missing(sigma)){
	if(is.vector(mean))  sigma <- as.vector(diag(length(mean)))
	else  sigma <- as.vector(diag(ncol(mean)))}
  if(is.vector(mean))  mean <- matrix(mean,nrow=1)
  nv <- ncol(mean) 
  if(nrow(mean) != n)  mean <- matrix(t(mean), ncol=nv, nrow=n, byrow=TRUE)

  if(is.vector(sigma)) {                       # 'classic' rmvnorm to gain time
    if(length(sigma) != (nv^2)) stop("sigma should be a vector of length:",nv^2)
    return(mean + rmvnorm(n, mean = rep(0, nv), sigma = matrix(sigma,ncol=nv), method=method))
        }

# else Sigma is a matrix

    if ((nv^2) != ncol(sigma)) stop("the length/ncol of sigma should be", nv^2, "\n")
    if (nrow(sigma) != n) sigma <- matrix(t(sigma), ncol = nv^2, nrow = n, byrow = TRUE)
    return(t(sapply(1:n, function(x) rmvnorm(1, mean[x, ], matrix(sigma[x,], ncol = nv), method = method))))
 }
 
#<<BEGIN>>
dmultinormal <- function (x, mean, sigma, log = FALSE) 
#ISALIAS rmultinormal
#--------------------------------------------

{
  # 'classic' dmvnorm to gain time
    if (is.vector(x))   x <- matrix(x, nrow = 1)
	nv <- ncol(x)
	if (missing(mean))  mean <- rep(0, length = nv)
    if (missing(sigma)) sigma <- as.vector(diag(nv))

	if(is.vector(mean) & is.vector(sigma)) return(dmvnorm(x=x, mean=mean, sigma = matrix(sigma, ncol=nv), log= log))

	if(is.vector(mean))  mean <- matrix(mean,nrow=1)
	if(is.vector(sigma)) sigma <- matrix(sigma,nrow=1)                      
    if(ncol(sigma) != (nv^2))  stop("the length/ncol of sigma should be", nv^2, "\n")
	
	if(nrow(mean) != nrow(x))  mean <- matrix(t(mean), ncol=nv, nrow=nrow(x), byrow=TRUE)
	if(nrow(sigma) != nrow(x)) sigma <- matrix(t(sigma), ncol=nv^2, nrow=nrow(x), byrow=TRUE)
	
    return(t(sapply(1:nrow(x), function(y) dmvnorm(x=x[y,], mean=mean[y, ], sigma=matrix(sigma[y,], ncol = nv), log = log))))
 }
