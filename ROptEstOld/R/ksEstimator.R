###############################################################################
## Kolmogorov(-Smirnov) minimum distance estimator
###############################################################################
setMethod("ksEstimator", signature(x = "numeric", distribution = "Binom"),
    function(x, distribution, param, eps = .Machine$double.eps^0.5){
        if(missing(param)){
            KSdist1 <- function(prob, x, size){
                supp <- 0:size
                edf <- ecdf(x)
                return(max(abs(edf(supp)-pbinom(supp, size = size, prob = prob))))
            }
            KSdist2 <- function(size, x, eps){
              unlist(optimize(f = KSdist1, interval = c(0, 1), tol = eps, 
                       x = x, size = size))
            }
            size <- max(x):(max(x) + 100)
            res <- sapply(size, KSdist2, x = x, eps = eps)
            ind <- which.min(res[2,])

            return(list(size = size[ind], prob = res[1,ind]))
        }
        if(param == "prob"){
            if(max(x) > size(distribution))
                stop("maximum of 'x' > 'size' of distribution")
            KSdist <- function(prob, x, size){
                supp <- 0:size
                edf <- ecdf(x)
                return(max(abs(edf(supp)-pbinom(supp, size = size, prob = prob))))
            }
            res <- optimize(f = KSdist, interval = c(0, 1), tol = eps, 
                     x = x, size = size(distribution))$minimum
            return(list(size = size(distribution), prob = res))
        }
        if(param == "size"){
            KSdist <- function(size, x, prob){
                supp <- 0:size
                edf <- ecdf(x)
                return(max(abs(edf(supp)-pbinom(supp, size = size, prob = prob))))
            }
            size <- max(x):(max(x) + 100)
            ind <- which.min(sapply(size, KSdist, x = x, prob = prob(distribution)))

            return(list(size = size[ind], prob = prob(distribution)))
        }
        stop("wrong 'param' specified")
    })

###############################################################################
## Minimum KS-distance for Poisson distribution
###############################################################################
setMethod("ksEstimator", signature(x = "numeric", distribution = "Pois"),
    function(x, distribution, param, eps = .Machine$double.eps^0.5){
        KSdist <- function(lambda, x){
            edf <- ecdf(x)
            supp <- 0:max(qpois(1-1e-8, lambda=lambda), max(x))
            return(max(abs(edf(supp)-ppois(supp, lambda=lambda))))
        }
        res <- optimize(f = KSdist, interval = c(0, max(x)), tol = eps, x = x)$minimum

        return(list(lambda = res))
    })

###############################################################################
## Minimum KS-distance for normal distribution
###############################################################################
setMethod("ksEstimator", signature(x = "numeric", distribution = "Norm"),
    function(x, distribution, param, eps = .Machine$double.eps^0.5){
        if(missing(param)){
            KSdist <- function(param, x){
                if(param[2] <= 0) return(Inf)
                return(ks.test(x, "pnorm", mean = param[1], sd = param[2])$statistic)
            }
            res <- optim(c(mean(distribution), sd(distribution)), f = KSdist, 
                        method = "Nelder-Mead", control=list(reltol=eps), 
                        x = x)$par
            return(list(mean = res[1], sd = res[2]))
        }
        if(param == "mean"){
            KSdist <- function(mean, x, sd){
                return(ks.test(x, "pnorm", mean = mean, sd = sd)$statistic)
            }
            res <- optimize(f = KSdist, interval = c(min(x), max(x)), 
                        tol = eps, x = x, sd = sd(distribution))$minimum
            return(list(mean = res, sd = sd(distribution)))
        }
        if(param == "sd"){
            KSdist <- function(sd, x, mean){
                return(ks.test(x, "pnorm", mean = mean, sd = sd)$statistic)
            }
            res <- optimize(f = KSdist, 
                        interval = c(.Machine$double.eps^0.5, max(x)-min(x)), 
                        tol = eps, x = x, mean = mean(distribution))$minimum
            return(list(mean = mean(distribution), sd = res))
        }
        stop("wrong 'param' specified")
    })

###############################################################################
## Minimum KS-distance for lognormal distribution
###############################################################################
setMethod("ksEstimator", signature(x = "numeric", distribution = "Lnorm"),
    function(x, distribution, param, eps = .Machine$double.eps^0.5){
        if(missing(param)){
            KSdist <- function(param, x){
                if(param[2] <= 0) return(Inf)
                return(ks.test(x, "plnorm", meanlog = param[1], sdlog = param[2])$statistic)
            }
            res <- optim(c(meanlog(distribution), sdlog(distribution)), f = KSdist, 
                        method = "Nelder-Mead", control=list(reltol=eps), 
                        x = x)$par
            return(list(meanlog = res[1], sdlog = res[2]))
        }
        if(param == "meanlog"){
            KSdist <- function(meanlog, x, sdlog){
                return(ks.test(x, "plnorm", meanlog = meanlog, sdlog = sdlog)$statistic)
            }
            res <- optimize(f = KSdist, interval = c(min(x), max(x)), 
                        tol = eps, x = x, sdlog = sdlog(distribution))$minimum
            return(list(meanlog = res, sdlog = sdlog(distribution)))
        }
        if(param == "sdlog"){
            KSdist <- function(sdlog, x, meanlog){
                return(ks.test(x, "plnorm", meanlog = meanlog, sdlog = sdlog)$statistic)
            }
            res <- optimize(f = KSdist, 
                        interval = c(.Machine$double.eps^0.5, max(x)-min(x)), 
                        tol = eps, x = x, meanlog = meanlog(distribution))$minimum
            return(list(meanlog = meanlog(distribution), sdlog = res))
        }
        stop("wrong 'param' specified")        
    })

###############################################################################
## Minimum KS-distance for Gumbel distribution
###############################################################################
setMethod("ksEstimator", signature(x = "numeric", distribution = "Gumbel"),
    function(x, distribution, param, eps = .Machine$double.eps^0.5){
        if(missing(param)){
            KSdist <- function(param, x){
                if(param[2] <= 0) return(Inf)
                return(ks.test(x, "pgumbel", loc = param[1], scale = param[2])$statistic)
            }
            res <- optim(c(mean(distribution), sd(distribution)), f = KSdist, 
                        method = "Nelder-Mead", control=list(reltol=eps), 
                        x = x)$par
            return(list(loc = res[1], scale = res[2]))
        }
        if(param == "loc"){
            KSdist <- function(loc, x, scale){
                return(ks.test(x, "pgumbel", loc = loc, scale = scale)$statistic)
            }
            res <- optimize(f = KSdist, interval = c(min(x), max(x)), 
                        tol = eps, x = x, scale = scale(distribution))$minimum
            return(list(loc = res, scale = scale(distribution)))
        }
        if(param == "scale"){
            KSdist <- function(scale, x, loc){
                return(ks.test(x, "pgumbel", loc = loc, scale = scale)$statistic)
            }
            res <- optimize(f = KSdist, 
                        interval = c(.Machine$double.eps^0.5, max(x)-min(x)), 
                        tol = eps, x = x, loc = loc(distribution))$minimum
            return(list(loc = loc(distribution), scale = res))
        }
        stop("wrong 'param' specified")
    })

###############################################################################
## Minimum KS-distance for Exponential distribution
###############################################################################
setMethod("ksEstimator", signature(x = "numeric", distribution = "Exp"),
    function(x, distribution, param, eps = .Machine$double.eps^0.5){
        KSdist <- function(rate, x){
            return(ks.test(x, "pexp", rate = rate)$statistic)
        }
        res <- optimize(f = KSdist, interval = c(0, max(x)), tol = eps, x = x)$minimum
        return(list(rate = res))
    })

###############################################################################
## Minimum KS-distance for Gamma model
###############################################################################
setMethod("ksEstimator", signature(x = "numeric", distribution = "Gammad"),
    function(x, distribution, param, eps = .Machine$double.eps^0.5){
        if(missing(param)){
            KSdist <- function(param, x){
                if((param[1] <= 0) || (param[2] <= 0)) return(Inf)
                return(ks.test(x, "pgamma", scale = param[1], shape = param[2])$statistic)
            }
            res <- optim(c(scale(distribution), shape(distribution)), f = KSdist, method = "Nelder-Mead", 
                        control=list(reltol=eps), x = x)$par
            return(list(scale = res[1], shape = res[2]))
        }
        if(param == "scale"){
            KSdist <- function(scale, x, shape){
                return(ks.test(x, "pgamma", scale = scale, shape = shape)$statistic)
            }
            res <- optimize(f = KSdist, interval = c(min(x), max(x)), 
                        tol = eps, x = x, shape = shape(distribution))$minimum
            return(list(scale = res, shape = shape(distribution)))
        }
        if(param == "shape"){
            KSdist <- function(shape, x, scale){
                return(ks.test(x, "pgamma", scale = scale, shape = shape)$statistic)
            }
            res <- optimize(f = KSdist, 
                        interval = c(.Machine$double.eps^0.5, max(x)-min(x)), 
                        tol = eps, x = x, scale = scale(distribution))$minimum
            return(list(scale = scale(distribution), shape = res))
        }
        stop("wrong 'param' specified")
    })
