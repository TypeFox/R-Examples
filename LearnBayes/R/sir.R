sir=function (logf, tpar, n, data) 
{
    k = length(tpar$m)
    theta = rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)
    lf=matrix(0,c(dim(theta)[1],1))
    for (j in 1:dim(theta)[1]) lf[j]=logf(theta[j,],data)
    lp = dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df, 
        log = TRUE)
    md = max(lf - lp)
    wt = exp(lf - lp - md)
    probs = wt/sum(wt)
    indices = sample(1:n, size = n, prob = probs, replace = TRUE)
    if (k > 1) 
        theta = theta[indices, ]
    else theta = theta[indices]
    return(theta)
}
