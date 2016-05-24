normchi2post=function(theta,data)
{
    mu = theta[1]
    sig2 = theta[2]
   
    logf=function(y,mu,sig2)
       -(y-mu)^2/2/sig2-log(sig2)/2

    z=sum(logf(data,mu,sig2))
    z = z - log(sig2)
    return(z)
}
