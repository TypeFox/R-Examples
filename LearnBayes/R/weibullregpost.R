weibullregpost=function (theta, data) 
{
    logf=function(t,c,x,sigma,mu,beta)
    {
    z=(log(t)-mu-x%*%beta)/sigma
    f=1/sigma*exp(z-exp(z))
    S=exp(-exp(z))
    c*log(f)+(1-c)*log(S)
    }
    
    k = dim(data)[2]
    p = k - 2
    t = data[, 1]
    c = data[, 2]
    X = data[, 3:k]
    sigma = exp(theta[1])
    mu = theta[2]
    beta = array(theta[3:k], c(p,1))
    return(sum(logf(t,c,X,sigma,mu,beta)))
}

