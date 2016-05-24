transplantpost=function (theta, data) 
{
    x = data[, 1]
    y = data[, 3]
    t = data[, 2]
    d = data[, 4]
    tau = exp(theta[1])
    lambda = exp(theta[2])
    p = exp(theta[3])

    xnt = x[t == 0]
    dnt = d[t == 0]
    z = x[t == 1]
    y = y[t == 1]
    dt = d[t == 1]

    logf=function(xnt,dnt,lambda,p)
      (dnt==0)*(p*log(lambda)+log(p)- (p + 1) * log(lambda + xnt)) +
      (dnt==1)*p*log(lambda/(lambda + xnt))

    logg=function(z,y,tau,lambda,p)
      (dt==0)*(p * log(lambda) + 
       log(p * tau)-(p + 1) * log(lambda + y + tau * z)) + 
        (dt==1) * p * log(lambda/(lambda + y + tau * z))

    val=sum(logf(xnt,dnt,lambda,p))+sum(logg(z,y,tau,lambda,p))
    
    val = val + theta[1] + theta[2] + theta[3]
    return(val)
}


