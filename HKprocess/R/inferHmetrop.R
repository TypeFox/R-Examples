inferHmetrop <- function(data,theta.init = 0.7,burnin = 500,mcmc = 20000,
                         thin = 1,tune = 1,verbose = 0,seed = NA) {
    # Set natural logarithm of eq. 10 in Tyralis and Koutsoyiannis (2014)
    logpHx <- function(H,x) {
        if (H <= 0) return(-Inf)
        if (H >= 1) return(-Inf)
        nx <- as.integer(length(x))
        maxlag <- nx - 1
        r <- acfHKp(H,maxlag)
        EPS <- .Machine$double.eps
        out <- .C("logpHx",as.double(r),nx,as.double(x),nx,EPS,
        tr = array(0,dim=c(1,5)),fault = as.integer(1),PACKAGE = "HKprocess")

        return(c((out$tr)[1],(out$tr)[2],(out$tr)[3],(out$tr)[4],(out$tr)[5]))
    }
    f <- function(H,x) (logpHx(H,x))[1]
    # Perform MCMC simulation
    
    return(MCMCmetrop1R(f,theta.init = theta.init,x = data,burnin = burnin,
    mcmc = mcmc,thin = thin,tune = tune,verbose = verbose,seed = seed))
}