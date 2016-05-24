# Description: 	A function to calculate posterior quantities of the bivariate normal.  See pages 79-86.
# Usage: 	biv.norm.post(data.mat,alpha,beta,m,n0=5) 
# Arguments: 	data.mat	A matrix with two columns of normally distributed data
#		alpha		Wishart first (scalar) parameter
#		beta 		Wishart second (matrix) parameter
#		m		prior mean for mu
#		n0		prior confidence parameter
# Values:	mu1		posterior mean, dimension 1
#		mu2		posterior mean, dimension 2
#		sig1		posterior variance, dimension 1
#		sig2		posterior variance, dimension 2
#		rho		posterior covariance



biv.norm.post <- function(data.mat,alpha,beta,m,n0=5) {
    n <- nrow(data.mat)
    xbar <- apply(data.mat,2,mean)
    S2 <- (n-1)*var(data.mat)
    Wp.inv <- solve(beta)+S2+((n0*n)/(n0+n))*(xbar-m)%*%t(xbar-m)
    Sigma <- solve(rwishart(alpha+n,SqrtSigma=chol(solve(Wp.inv))))
    mu <- rmultinorm(1, (n0*m + n*xbar)/(n0+n), Sigma/(n0+n))
    return(c(mu1=mu[1],mu2=mu[2],sig1=Sigma[1,1], sig2=Sigma[2,2],rho=Sigma[2,1]))
}






