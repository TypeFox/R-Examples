betabinexch0=function (theta, data) 
{
    eta = theta[1]
    K = theta[2]
    y = data[, 1]
    n = data[, 2]
    N = length(y)

    logf=function(y,n,K,eta)
       lbeta(K * eta + y, K * (1 - eta) + n - y)-lbeta(K * eta, K * (1 - eta))

    val=sum(logf(y,n,K,eta))

    val = val - 2 * log(1 + K) - log(eta) - log(1 - eta)
    return(val)
}
