pmmglm <- function(x, beta, sigma, glmfamily, Xdesign, size=NA, log.p=FALSE){
    mu <- glmfamily$linkinv(Xdesign %*% beta)
    if (glmfamily$family=="gaussian")
        return(pnorm(x, mean=mu, sd=sigma, log.p=log.p))
    else if (glmfamily$family=="poisson")
        return(ppois(x, lambda=mu, log.p=log.p))
    else if (glmfamily$family=="Gamma")
        return(pgamma(x, scale=mu*sigma^2, shape=1/sigma^2, log.p=log.p))
    else if (glmfamily$family=="binomial")
        return(pbinom(x, size=size, prob=mu, log.p=log.p))
}


