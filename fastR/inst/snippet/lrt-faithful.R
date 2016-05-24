# density function for mixture of normals
dmix <- function(x, alpha,mu1,mu2,sigma1,sigma2) {
    if (alpha < 0) return (dnorm(x,mu2,sigma2))
    if (alpha > 1) return (dnorm(x,mu1,sigma1))
    
    alpha * dnorm(x,mu1,sigma1) + (1-alpha) * dnorm(x,mu2,sigma2)
    }

# log-likelihood
loglik <- function(theta, x) {
    alpha <- theta[1]  
    mu1 <- theta[2]
    mu2 <- theta[3]
    sigma1 <- theta[4]
    sigma2 <- theta[5]
    density <- function (x) {
        if (alpha < 0) return (Inf)
        if (alpha > 1) return (Inf)
        if (sigma1<= 0) return (Inf)
        if (sigma2<= 0) return (Inf)
        dmix(x,alpha,mu1,mu2,sigma1,sigma2)
    }
    sum( log ( sapply( x, density) ) )
}

loglik0 <- function(theta, x) {
    theta <- c(0.5,theta)
    return(loglik(theta,x))
}
# seed the algorithm  
m <- mean(faithful$eruptions)
s <- sd(faithful$eruptions)

oldopt <- options(warn=-1)    # suppress warnings from log(0) 
mle <-  nlmax(loglik, p=c(0.5,m-1,m+1,s,s), x=faithful$eruptions)$estimate
mle
loglik(mle,x=faithful$eruptions)
mle0 <- nlmax(loglik0,p=c(m-1,m+1,s,s), x=faithful$eruptions)$estimate
mle0
loglik0(mle0,x=faithful$eruptions)
stat <- 2 * (loglik(mle,x=faithful$eruptions) 
           - loglik0(mle0,x=faithful$eruptions)); stat
1 - pchisq(stat,df=1)         # p-value based on asymptotic distribution
options(oldopt)
