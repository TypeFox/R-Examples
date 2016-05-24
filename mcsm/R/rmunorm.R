rmunorm=function(mu,sig){
#generator from the multivariate Gaussian

as.vector(mu+t(chol(sig))%*%rnorm(length(mu)))
}
