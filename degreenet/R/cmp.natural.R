cmp.mutonatural <- function (mu,sig,K=1000) 
{
#   natural is log(lambda), -nu
    options(warn = -1)
    result = optim(par=c(log(mu), 0), fn=function(p,mu,sig,K) {
        j <- 0:K
        a <- dcmp.natural(x=j, v=c(exp(p[1]), exp(p[2])), err=0.00000001)
	sqrt(abs(sum(a*j)^2-mu*mu)+abs(sum(a*j*j)-sig*sig-mu*mu))
    }, mu=mu, sig=sig, K=K, method = "Nelder-Mead",
       control=list(maxit=20000))
    options(warn = 0)
    lambda = exp(result$par[1])
    nu = exp(result$par[2])
    fit = c(lambda = lambda, nu = nu, result)
    return(fit)
}
cmp.naturaltomu <- function (p,K=1000) 
{
        j <- 0:K
        a <- dcmp.natural(x=j, v=p, err=0.00000001)
	mu <- sum(a*j)
	sd <- sqrt(sum(a*j*j)-mu*mu)
    return(c(mu,sd))
}
