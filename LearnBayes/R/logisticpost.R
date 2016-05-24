logisticpost=function (beta, data) 
{
    x = data[, 1]
    n = data[, 2]
    y = data[, 3]

    beta0 = beta[1]
    beta1 = beta[2]
    
    logf=function(x,n,y,beta0,beta1)
    {  lp = beta0 + beta1 * x
       p = exp(lp)/(1 + exp(lp))
       y * log(p) + (n - y) * log(1 - p)
     }

    return(sum(logf(x,n,y,beta0,beta1)))
}
