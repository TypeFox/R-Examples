oldopt <- options(warn=-1)          # suppress warnings from log(0) 
theta2probs <- function(theta) { 
    c(theta^2, 2*theta*(1-theta), (1-theta)^2)  
}
loglik <- function(theta,x) {
    probs <- theta2probs(theta)
    if (any (probs <=0)) return (Inf)
    return ( dmultinom(x,size=sum(x),prob=probs,log=T) )
}
geno<-c(83,447,470)
nlmax(loglik,p=0.5,x=geno)$estimate
options(oldopt)                         # reset options
