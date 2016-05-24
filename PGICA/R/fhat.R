fhat <-
function(x,theta,mu,sigma){sum(theta*dnorm(x,mu,sigma))}
