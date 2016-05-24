dpgnorm <-
function(y,p,mean,sigma)
{
# A function implemented by Steve Kalke

# Description: 
# Computes the density function of the p-generalized normal distribution 
# for the real Argument "y"

# Arguments: 
# p- a positiv constant (default: p=2)
# mean- a real constant, expressing the expectation (default: mean=0)
# sigma- a positiv constant, expressing the standard deviation (default: sigma=1)

if(missing(mean)){mean<-0}

if(missing(sigma)){sigma<-1}

if(missing(p)){p<-2}

if(p<0){stop("p has to be positive")}

x<-(y-mean)/sigma

density<- 1/(2*p^(1/p-1)*gamma(1/p))*exp(-(abs(x)^p)/p)

return(density)

}
