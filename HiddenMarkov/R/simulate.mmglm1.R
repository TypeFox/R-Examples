simulate.mmglm1 <- function (object, nsim=1, seed=NULL, ...){
    if (nsim!=1)
        warning("Argument 'nsim' is redundant, determined by nrow(Xdesign)")
    nsim <- nrow(object$Xdesign)
    if (!is.null(seed)) set.seed(seed)
    mstate <- simulate.mchain(object, nsim=nsim, seed=seed)$mc
    linkinv <- object$glmfamily$linkinv
    distribution <- object$glmfamily$family
    #   redundancy below, may be faster than loop though
    mu <- linkinv(object$Xdesign %*% object$beta)
    mu <- mu[cbind(1:nsim, mstate)]
    if (distribution=="gaussian")
        object$y <- rnorm(nsim, mean=mu, sd=object$sigma[mstate])
    else if (distribution=="poisson")
        object$y <- rpois(nsim, lambda=mu)
    else if (distribution=="Gamma")
        object$y <- rgamma(nsim, scale=mu*object$sigma[mstate]^2,
                           shape=1/object$sigma[mstate]^2)
    else if (distribution=="binomial")
        object$y <- rbinom(nsim, size=object$size, prob=mu)
    return(object)
}


