truncnorm=function(n,mu,tau2,a,b)
{
qnorm(pnorm(b,mu,sqrt(tau2))-runif(n)*(pnorm(b,mu,sqrt(tau2))-pnorm(a,mu,sqrt(tau2))),mu,sqrt(tau2))
}
