pglm <- function(q, x1, beta0, beta1, sigma, family, link, size=NA, log.p=FALSE){
    eta <- beta0 + beta1*x1
    #-----------------------------
    if (link=="inverse") mu <- 1/eta
    else if (link=="identity") mu <- eta
    else if (link=="log") mu <- exp(eta)
    else if (link=="logit") prob <- exp(eta)/(1+exp(eta))
    else if (link=="probit") prob <- pnorm(eta)
    else if (link=="cloglog") prob <- 1-exp(-exp(eta))
    #-----------------------------
    if (family=="gaussian")
        return(pnorm(q, mean=mu, sd=sigma, log.p=log.p))
    else if (family=="poisson")
        return(ppois(q, lambda=mu, log.p=log.p))
    else if (family=="Gamma")
        return(pgamma(q, scale=mu*sigma^2, shape=1/sigma^2, log.p=log.p))
    else if (family=="binomial")
        return(pbinom(q, size=size, prob=prob, log.p=log.p))
}

