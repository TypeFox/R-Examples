###
### momunknown.R
###

momunknown <- function(theta1hat,V1,n,nuisance.theta,g=1,theta0,ssr,logbf=FALSE) {
if (missing(theta0)) theta0 <- rep(0,length(theta1hat))
p1 <- length(theta1hat); p <- p1 + nuisance.theta
l <- theta1hat-theta0; l <- matrix(l,nrow=1) %*% solve(V1) %*% matrix(l,ncol=1)  
sigma2hat <- (ssr + l/(1+n*g))/(n-nuisance.theta)
muk <- p1+ l* n*g/((1+n*g)*sigma2hat)
bf <- (-(n-nuisance.theta)/2)*log(1+n*g*ssr/(ssr+l)) + log(muk) - log(1+n*g) + ((n-p)/2)*log(1+n*g) - log(p1)
if (!logbf) bf <- exp(bf)
return(bf)
}

