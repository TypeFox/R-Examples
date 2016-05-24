normpostsim=function (data, prior=NULL, m = 1000) 
{

if (length(prior)==0)
{
    S = sum((data - mean(data))^2)
    xbar = mean(data)
    n = length(data)
    SIGMA2 = S/rchisq(m, n - 1)
    MU = rnorm(m, mean = xbar, sd = sqrt(SIGMA2)/sqrt(n))
} else
{
    a=prior$sigma2[1]
    b=prior$sigma2[2]
    mu0=prior$mu[1]
    tau2=prior$mu[2]
    S = sum((data - mean(data))^2)
    xbar = mean(data)
    n = length(data)
 
    SIGMA2=rep(0,m)
    MU=rep(0,m)
    sigma2=S/n
    for (j in 1:m)
    {
    prec=n/sigma2+1/tau2
    mu1=(xbar*n/sigma2+mu0/tau2)/prec
    v1=1/prec
    mu=rnorm(1,mu1,sqrt(v1))

    a1=a+n/2
    b1=b+sum((data-mu)^2)/2
    sigma2=rigamma(1,a1,b1)

    SIGMA2[j]=sigma2
    MU[j]=mu
    }
}
    return(list(mu = MU, sigma2 = SIGMA2))
}
