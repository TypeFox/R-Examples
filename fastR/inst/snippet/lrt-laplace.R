x <- c(1.00,-1.43,0.62,0.87,-0.66,-0.59,1.30,-1.23,-1.53,-1.94)
loglik1 <- function(theta, x) {
    m <- theta[1]; lambda <- theta[2]
    return( sum( log(0.5) + dexp(abs(x-m),rate=lambda, log=T)) )
}
loglik0 <- function(theta, x) {
    m <- 0; lambda <- theta[1]
    return( sum( log(0.5) + dexp(abs(x-m),rate=lambda, log=T)) )
}
oldopt <- options(warn=-1)
free <- nlmax(loglik1,p=c(0,1),x=x)$estimate; free
null <- nlmax(loglik0,p=c(1),x=x)$estimate; null
stat <- 2 * (loglik1(free,x) - loglik0(null,x)); stat
1 - pchisq(stat,df=1)          # p-value based on asymptotic distribution
options(oldopt)
