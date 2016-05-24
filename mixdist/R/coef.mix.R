## last modified May 2008

coef.mix <- function(object, natpar = FALSE, ...) 
{
    mixobj<-object
    par <- mixobj$parameters
    dist <- mixobj$distribution
    constr <- mixobj$constraint
    mu <- par[, 2]
    sigma <- par[, 3]
    if (!natpar) 
        coef <- par
    else {
        if (dist == "norm") 
            coef <- par
        else if (dist == "lnorm") {
            scale <- sqrt(log((sigma/mu)^2 + 1))
            shape <- log(mu) - (scale^2)/2
            coef <- cbind(par, shape, scale)
        }
        else if (dist == "gamma") {
            shape <- (mu/sigma)^2
            rate <- mu/(sigma^2)
            coef <- cbind(par, shape, rate)
        }
        else if (dist == "weibull") {
            weibpar <- weibullpar(mu, sigma)
            shape <- weibpar$shape
            scale <- weibpar$scale
            coef <- cbind(par, shape, scale)
        }
        else if (dist == "binom") {
            size <- constr$size
            prob <- mu/constr$size
            coef <- cbind(par, size, prob)
        }
        else if (dist == "nbinom") {
            if (constr$consigma == "NBINOM") 
                size <- constr$size
            else size <- mu^2/(sigma^2 - mu)
            prob <- size/(size + mu)
            coef <- cbind(par, size, prob)
        }
        else if (dist == "pois") {
            lambda <- mu
            coef <- cbind(par, lambda)
        }
    }
    coef
}
