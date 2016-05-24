oldopt <- options(warn=-1)    # suppress warnings from log(0) 
# loglik defined above
mle <-  nlmax(loglik, p=c(0.5,m-1,m+1,s,s), x=faithful$eruptions)$estimate
f <- function(a) {
    loglik0 <- function(theta, x) {
        theta <- c(a,theta)
        return(loglik(theta,x))              
    }
    mle0 <- nlmax(loglik0,p=c(m-1,m+1,s,s), x=faithful$eruptions)$estimate
    stat <- 2 * (loglik(mle,x=faithful$eruptions) 
               - loglik0(mle0,x=faithful$eruptions)); stat
    pval <- 1 - pchisq(stat,df=1)         
	return(pval)
}
uniroot( function(a){f(a) - 0.05}, c(0,mle[1]))$root
uniroot( function(a){f(a) - 0.05}, c(1,mle[1]))$root
options(oldopt)
