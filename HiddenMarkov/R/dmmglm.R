dmmglm <- function(x, beta, sigma, glmfamily, Xdesign, size=NA, log=FALSE){
    mu <- glmfamily$linkinv(Xdesign %*% beta)
    if (glmfamily$family=="gaussian")
        return(dnorm(x, mean=mu, sd=sigma, log=log))
    else if (glmfamily$family=="poisson")
        return(dpois(x, lambda=mu, log=log))
    else if (glmfamily$family=="Gamma")
        return(dgamma(x, scale=mu*sigma^2, shape=1/sigma^2, log=log))
    else if (glmfamily$family=="binomial")
        return(dbinom(x, size=size, prob=mu, log=log))
}


