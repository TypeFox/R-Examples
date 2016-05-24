gaussprod<-function(mu1,mu2,sig1=1,sig2=1)
{
resp<-(2*pi)^(-1/2)*(sig1^2+sig2^2)^(-1/2)*exp(-(mu1-mu2)^2/(2*(sig1^2+sig2^2)))
return(resp)
}







