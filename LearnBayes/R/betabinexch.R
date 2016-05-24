betabinexch=function (theta, data) 
{
    eta = exp(theta[1])/(1 + exp(theta[1]))
    K = exp(theta[2])
    y = data[, 1]
    n = data[, 2]
    N = length(y)
 
    logf=function(y,n,K,eta)
       lbeta(K * eta + y, K * (1 - eta) + n - y)-lbeta(K * eta, K * (1 - eta))

    val=sum(logf(y,n,K,eta))

    val = val + theta[2] - 2 * log(1 + exp(theta[2]))
    return(val)
}
