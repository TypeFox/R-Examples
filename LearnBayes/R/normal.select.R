normal.select=function (quantile1, quantile2) 
{
    p1 = quantile1$p
    x1 = quantile1$x
    p2 = quantile2$p
    x2 = quantile2$x
    
    sigma=(x1-x2)/diff(qnorm(c(p2,p1)))
    mu=x1-sigma*qnorm(p1)

    return(list(mu=mu,sigma=sigma))
}
