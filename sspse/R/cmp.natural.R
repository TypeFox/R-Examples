cmp.natural <- function (mu,sigma,K=max(ceiling(mu+10*sigma),10),guess=NULL) 
{
# converts mean to natural
# natural is log(lambda), -nu
# natural is lambda, nu
#   K <- max(ceiling(mu+10*sigma),10)
    options(warn = -1)
    if(is.null(guess)){
      guess = c(2*log(mu+0.25), log(2))
    }
    result = optim(par=guess, fn=function(p,mu,sigma,K) {
        j <- 0:K
        a <- dcmp(x=j, lambda=exp(p[1]), nu=exp(p[2]), err=0.000000001)
	sqrt(abs(sum(a*j)^2-mu*mu)+abs(sum(a*j*j)-sigma*sigma-mu*mu))
    }, mu=mu, sigma=sigma, K=K, method = "Nelder-Mead",
       control=list(maxit=20000))
    options(warn = 0)
    lambda = exp(result$par[1])
    nu = exp(result$par[2])
    fit = c(lambda = lambda, nu = nu, result)
    return(fit)
}
cmp.mu <- function (p,K=1000,cutoff=1,max.mu=50) 
{
# converts natural to mean
        j <- cutoff:K
#       a <- dcmp(x=j, lambda=p[1], nu=p[2], err=0.000000001)
        a <- j*log(p[1])-p[2]*lgamma(j+1.0)
        a <- a - max(a,na.rm=TRUE)
        a[is.na(a)] <- min(a,na.rm=TRUE)
        a <- exp(a)
        a <- a / sum(a)
	mu <- sum(a*j)
	sd <- sqrt(sum(a*j*j)-mu*mu)
#   cat(sprintf("K= %d mu= %f sd= %f\n",K,mu,sd))
        K <- max(ceiling(mu+10*sd),10)
        j <- cutoff:K
#       a <- dcmp(x=j, lambda=p[1], nu=p[2], err=0.000000001)
        a <- j*log(p[1])-p[2]*lgamma(j+1.0)
        a <- a - max(a,na.rm=TRUE)
        a[is.na(a)] <- min(a,na.rm=TRUE)
        a <- exp(a)
        a <- a / sum(a)
	mu <- sum(a*j)
	sd <- sqrt(sum(a*j*j)-mu*mu)
#   cat(sprintf("K= %d mu= %f sd= %f\n",K,mu,sd))
        if(mu > max.mu){
          return(c(NA,NA))
        }else{
          return(c(mu,sd))
        }
}
