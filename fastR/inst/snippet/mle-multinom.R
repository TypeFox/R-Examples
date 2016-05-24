oldopt <- options(warn=-1)    # suppress warnings from log(0) 
loglik <- function(theta,x) { 
    probs <- c(theta, 1-sum(theta))
    if (any (probs < 0)) {return(Inf)}
    return( dmultinom(x,size=100,prob=probs, log=T) )
    }
nlmax(loglik,p=rep(0.25,3),x=c(10,20,30,40))$estimate -> mle; mle
options(oldopt)                     # restore options
round(mle,6)
