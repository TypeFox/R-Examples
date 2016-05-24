dmunorm=function(x,mu,sig,log=FALSE){
#density of the multivariate Gaussian

vale=-.5*((x-mu)%*%solve(sig)%*%(x-mu)+log(det(sig))+length(x)*log(2*pi))
if (!log) vale=exp(vale)

vale
}
