simulate.mmglm0 <- function (object, nsim=1, seed=NULL, ...){
    if (!is.null(seed)) set.seed(seed)
    obj1 <- simulate.mchain(object, nsim=nsim, seed=seed)
    #   simulate covariate as uniform (0,1) if NULL
    if (is.null(obj1$x$x1)) obj1$x$x1 <- runif(nsim, min=0, max=1)
    else if (length(obj1$x$x1)!=nsim) stop("length(x1) ne nsim")
    #   used for binomial case: number of Bernoulli trials
    if (obj1$family=="binomial"){
        if (is.null(obj1$x$size)) obj1$x$size <- 100+rpois(nsim, lambda=5)
        else if (length(obj1$x$size)!=nsim) stop("length(xsize) ne nsim")
    }
    obj1$x$y <- rep(NA,nsim)
    for (i in 1:nsim){
        eta <- obj1$beta[1,obj1$mc[i]] + obj1$beta[2,obj1$mc[i]]*obj1$x$x1[i]
        #---------------------------------------------
        if (obj1$link=="inverse") mu <- 1/eta
        else if (obj1$link=="identity") mu <- eta
        else if (obj1$link=="log") mu <- exp(eta)
        else if (obj1$link=="logit") prob <- exp(eta)/(1+exp(eta))
        else if (obj1$link=="probit") prob <- pnorm(eta)
        else if (obj1$link=="cloglog") prob <- 1-exp(-exp(eta))
        #---------------------------------------------
        if (obj1$family=="gaussian")
            obj1$x$y[i] <- rnorm(1, mean=mu, sd=obj1$sigma[obj1$mc[i]])
        else if (obj1$family=="poisson")
            obj1$x$y[i] <- rpois(1, lambda=mu)
        else if (obj1$family=="Gamma")
            obj1$x$y[i] <- rgamma(1, scale=mu*obj1$sigma[obj1$mc[i]]^2,
                             shape=1/obj1$sigma[obj1$mc[i]]^2)
        else if (obj1$family=="binomial")
            obj1$x$y[i] <- rbinom(1, size=obj1$x$size[i], prob=prob)
    }
    obj1$x <- as.data.frame(obj1$x)
    return(obj1)
}

