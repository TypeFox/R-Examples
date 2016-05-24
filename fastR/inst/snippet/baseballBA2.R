# log likelihood function
loglik <- function(theta,x) { sum(dbeta(x,theta[1],theta[2],log=T)) }
