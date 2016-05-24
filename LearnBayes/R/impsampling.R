impsampling=function (logf, tpar, h, n, data) 
{
    theta = rmt(n, mean = c(tpar$m), S = tpar$var, df = tpar$df)

    lf=matrix(0,c(dim(theta)[1],1))
    for (j in 1:dim(theta)[1]) lf[j]=logf(theta[j,],data)
    H=lf
    for (j in 1:dim(theta)[1]) H[j]=h(theta[j,])

    lp = dmt(theta, mean = c(tpar$m), S = tpar$var, df = tpar$df, 
        log = TRUE)
    md = max(lf - lp)
    wt = exp(lf - lp - md)
    est = sum(wt * H)/sum(wt)
    SEest = sqrt(sum((H - est)^2 * wt^2))/sum(wt)
    return(list(est = est, se = SEest, theta = theta, wt = wt))
}
