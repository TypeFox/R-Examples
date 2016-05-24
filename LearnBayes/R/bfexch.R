bfexch=function (theta, datapar) 
{
    y = datapar$data[, 1]
    n = datapar$data[, 2]
    K = datapar$K
    eta = exp(theta)/(1 + exp(theta))
   
    logf=function(K,eta,y,n)
         lbeta(K*eta+y, K*(1-eta)+n-y)-lbeta(K*eta, K*(1-eta))

    sum(logf(K,eta,y,n)) + log(eta * (1 - eta))-
         lbeta(sum(y) + 1, sum(n - y) + 1)
}
