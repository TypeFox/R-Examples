dglm <- function(x, x1, beta0, beta1, sigma, family, link, size=NA, log=FALSE){
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
        return(dnorm(x, mean=mu, sd=sigma, log=log))
    else if (family=="poisson")
        return(dpois(x, lambda=mu, log=log))
    else if (family=="Gamma")
        return(dgamma(x, scale=mu*sigma^2, shape=1/sigma^2, log=log))
    else if (family=="binomial")
        return(dbinom(x, size=size, prob=prob, log=log))
}

