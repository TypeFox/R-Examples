 poissgamexch=function (theta, datapar)
{
    y = datapar$data[, 2]
    e = datapar$data[, 1]
    z0 = datapar$z0
    alpha = exp(theta[1])
    mu = exp(theta[2])
    beta = alpha/mu
 
    logf=function(y,e,alpha,beta)
       lgamma(alpha + y) - (y + alpha) * log(e + beta) +
            alpha * log(beta)-lgamma(alpha)

    val=sum(logf(y,e,alpha,beta))
    val = val + log(alpha) - 2 * log(alpha + z0)
    return(val)
}