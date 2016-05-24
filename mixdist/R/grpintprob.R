## last modified June 2002

grpintprob <- function(mixdat, mixpar, dist, constr) 
{
    m <- nrow(mixdat)
    k <- nrow(mixpar)
    mu <- mixpar[, 2]
    sigma <- mixpar[, 3]
    if (dist == "norm") {
        par1 <- mu
        par2 <- sigma
        mixcdf <- t(sapply(mixdat[-m, 1], pnorm, par1, par2))
    }
    else if (dist == "lnorm") {
        par2 <- sqrt(log((sigma/mu)^2 + 1))
        par1 <- log(mu) - (par2^2)/2
        mixcdf <- t(sapply(mixdat[-m, 1], plnorm, par1, par2))
    }
    else if (dist == "gamma") {
        par1 <- (mu/sigma)^2
        par2 <- mu/(sigma^2)
        mixcdf <- t(sapply(mixdat[-m, 1], pgamma, par1, par2))
    }
    else if (dist == "weibull") {
        par <- weibullpar(mu, sigma)
        par1 <- par$shape
        par2 <- par$scale
        mixcdf <- t(sapply(mixdat[-m, 1], pweibull, par1, par2))
    }
    else if (dist == "binom") {
        par1 <- constr$size
        par2 <- mu/constr$size
        mixcdf <- t(sapply(mixdat[-m, 1], pbinom, par1, par2))
    }
    else if (dist == "nbinom") {
        if (constr$consigma == "NBINOM") 
            par1 <- constr$size
        else par1 <- mu^2/(sigma^2 - mu)
        mixcdf <- t(sapply(mixdat[-m, 1], pnbinom, par1, mu = mu))
    }
    else if (dist == "pois") {
        par <- mu
        mixcdf <- t(sapply(mixdat[-m, 1], ppois, par))
    }
    if (k == 1) 
        mixcdf <- t(mixcdf)
    rbind(mixcdf, 1) - rbind(0, mixcdf)
}
