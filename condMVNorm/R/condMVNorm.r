condMVN <- function(mean, sigma, dependent.ind, given.ind, X.given, check.sigma=TRUE){
# Returns conditional mean and variance of 
# X[dependent.ind] | X[given.ind] = X.given
# where X is multivariateNormal(mean = mean, covariance = sigma)
# 
#
if (missing(dependent.ind)) return("You must specify the indices of dependent random variables in `dependent.ind'")
if (missing(given.ind) & missing(X.given)) return(list(condMean=mean[dependent.ind], condVar=as.matrix(sigma[dependent.ind, dependent.ind])))
if (length(X.given) != length(given.ind)) stop("lengths of `X.given' and `given.ind' must be same")
if (check.sigma) {
	if (!isSymmetric(sigma)) 
        stop("sigma is not a symmetric matrix")
	eigenvalues <- eigen(sigma, only.values = TRUE)$values
	if (any(eigenvalues < 1.e-08))
        stop("sigma is not positive-definite")
}
B <- sigma[dependent.ind, dependent.ind]
C <- sigma[dependent.ind, given.ind, drop=FALSE]
D <- sigma[given.ind, given.ind]
CDinv <- C %*% solve(D)
cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
cVar <- B - CDinv %*% t(C)
list(condMean=cMu, condVar=cVar)
}

# PDF of conditional multivariate-normal
dcmvnorm <- function(x, mean, sigma, dependent.ind, given.ind, X.given, check.sigma=TRUE, log=FALSE){
if (missing(dependent.ind)) return("You must specify the indices of dependent random variables in `dependent.ind'")
if (length(x) != length(dependent.ind)) stop("lengths of `x' and `dependent.ind' must be same")
if (missing(given.ind) & missing(X.given)) return(dmvnorm(x, mean=mean[dependent.ind], sigma=as.matrix(sigma[dependent.ind, dependent.ind]), log=log))
if (length(X.given) != length(given.ind)) stop("lengths of `X.given' and `given.ind' must be same")
ret <- condMVN(X.given=X.given, mean=mean, sigma=sigma, dependent.ind=dependent.ind, given.ind=given.ind, check.sigma=check.sigma)
dmvnorm(x, mean=ret$condMean, sigma=ret$condVar, log=log)
}

# CDF of conditional multivariate-normal
pcmvnorm <- function(lower=-Inf, upper=Inf, mean, sigma, dependent.ind, given.ind, X.given, check.sigma=TRUE, algorithm=GenzBretz(), ...){
if (missing(dependent.ind)) return("You must specify the indices of dependent random variables in `dependent.ind'")
if (missing(given.ind) & missing(X.given)) return(pmvnorm(lower=lower, upper=upper, mean=mean[dependent.ind], sigma=as.matrix(sigma[dependent.ind, dependent.ind])))
if (length(X.given) != length(given.ind)) stop("lengths of `X.given' and `given.ind' must be same")
ret <- condMVN(X.given=X.given, mean=mean, sigma=sigma, dependent.ind=dependent.ind, given.ind=given.ind, check.sigma=check.sigma)
pmvnorm(lower=lower, upper=upper, mean=ret$condMean, sigma=ret$condVar, algorithm=algorithm, ...)
}

# Random-number generation from conditional multivariate-normal
rcmvnorm <- function(n, mean, sigma, dependent.ind, given.ind, X.given, check.sigma=TRUE, method=c("eigen", "svd", "chol")){
if (missing(dependent.ind)) return("You must specify the indices of dependent random variables in `dependent.ind'")
if (missing(given.ind) & missing(X.given)) return(rmvnorm(n,  mean=mean[dependent.ind], sigma=as.matrix(sigma[dependent.ind, dependent.ind])))
if (length(X.given) != length(given.ind)) stop("lengths of `X.given' and `given.ind' must be same")
ret <- condMVN(X.given=X.given, mean=mean, sigma=sigma, dependent.ind=dependent.ind, given.ind=given.ind, check.sigma=check.sigma)
rmvnorm(n, mean=ret$condMean, sigma=ret$condVar ,
        method=method)
}
