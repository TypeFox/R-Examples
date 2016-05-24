groupeddatapost=function (theta, data) 
{
    dj=function(f,int.lo,int.hi,mu,sigma)
      f*log(pnorm(int.hi,mu,sigma)-pnorm(int.lo,mu,sigma))
    mu = theta[1]
    sigma = exp(theta[2])
    sum(dj(data$f,data$int.lo,data$int.hi,mu,sigma))
}
