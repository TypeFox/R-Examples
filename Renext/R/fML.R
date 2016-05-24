##==============================================================================
## Author: Yves Deville
##
## Find ML estimates for special distributions
##
## All parameters (e.g. dist. name) must here be given in clean form.
##
##==============================================================================

fML <- function(y,
                distname.y,
                parnames.y,
                fixed.y,
                fixed.par.y) {
    
    ## cat("fixed.y\n")
    ## print(fixed.y)
    ## cat("fixed.par.y\n")
    ## print(fixed.par.y)
    
    n <- length(y)
    
    if(distname.y == "exponential") {
        
        ##======================================================================
        ## Explicit maximum likelihood
        ##======================================================================
        
        if (fixed.y) {
            rate.hat <- fixed.par.y
            cov0.y <- matrix(0, nrow = 1, ncol = 1)
            logLik <- n * log(rate.hat) - rate.hat * sum(y)
        } else {
            rate.hat <- 1 / mean(y)
            est.y <- rate.hat
            cov0.y <- matrix(rate.hat^2 / n, nrow = 1, ncol = 1)
            logLik <- n * (log(rate.hat) - 1.0)
        }
        
        names(est.y) <- parnames.y
        colnames(cov0.y) <- rownames(cov0.y) <- "rate"
        
    } else if(distname.y == "weibull") {
        
        ##======================================================================
        ## Explicit maximum likelihood for known shape,
        ## ML with concentration with 'fweibull' for the general case
        ##======================================================================
        
        if (fixed.y["shape"]) {
            
            shape.hat <- fixed.par.y["shape"]
            if (shape.hat <= 0) stop("'fixed.par.y' must specify a 'shape' > 0")
            
            if (fixed.y["scale"]) {
                scale.hat <- fixed.par.y["scale"]
                ymod <- y / scale.hat
                cov0.y <- matrix(0, nrow = 2, ncol = 2)
                logLik <- n * ( log(shape.hat / scale.hat) +
                               (shape.hat - 1) * mean(log(ymod)) -
                               mean(ymod^shape.hat) )
            } else {
                scale.hat <- mean(y^shape.hat)^{1/shape.hat}
                ymod <- y / scale.hat
                cov0.y <- matrix(0, nrow = 2, ncol = 2)
                cov0.y[2, 2] <- ( (scale.hat/shape.hat)^2 ) / n
                logLik <- n * ( log(shape.hat / scale.hat) +
                               (shape.hat - 1) * mean(log(ymod)) - 1.0)
            }
            
        } else if (fixed.y["scale"]) {
            stop("parameter 'scale' can not be fixed alone in Weibull")
        } else {
            ## concentrated maximum likelihood
            fit <- fweibull(x = y, info.observed = FALSE)
            est <- fit$estimate
            shape.hat <- est["shape"]
            scale.hat <- est["scale"]
            cov0.y <- fit$cov
            logLik <- fit$loglik
        }
        
        ## CAUTION: here respect the order specified above for Weibull!
        est.y <- c(shape.hat, scale.hat)
        names(est.y) <- parnames.y
        colnames(cov0.y) <- rownames(cov0.y) <- parnames.y
        
    } else if (distname.y =="log-normal") {
        
        ##=====================================================================
        ## Explicit ML
        ##=====================================================================
        
        ly <- log(y)
        meanlog.hat <- mean(ly)
        sdlog.hat <- sd(ly) * sqrt((n-1)/n)
    
        est.y <- c(meanlog.hat, sdlog.hat)
        names(est.y) <- parnames.y
        
        sig2 <- sdlog.hat^2
        ## Do not forget that 
        cov0.y <- matrix(c(sig2/n, 0, 0, sig2/2/n), nrow = 2L, ncol = 2L)
        colnames(cov0.y) <- rownames(cov0.y) <- names(est.y)
        
        logLik <- n* (-log(sdlog.hat * sqrt((2*pi))) -
                      meanlog.hat - sdlog.hat^2 )
        
    } else if(distname.y == "gpd") {
        
        ##=====================================================================
        ## Specific maximum likelihood for known shape (concave loglik),
        ## ML with 'evd::fpot' for the general case
        ##=====================================================================
        
        if (fixed.y["scale"]) {
            stop("parameter 'scale' can not be fixed alone in \"gpd\"")
        } else if (fixed.y["shape"]) {
            shape.hat <- fixed.par.y["shape"]
            fit <- fgpd1(x = y, shape = fixed.par.y["shape"])
            scale.hat <- fit$estimate
            est.y <- c(scale.hat, shape.hat)
            names(est.y) <- c("scale", "shape")
            cov0.y <- matrix(c(fit$cov, 0, 0, 0), nrow = 2L, ncol = 2L)
            colnames(cov0.y) <- rownames(cov0.y) <- names(est.y)
            logLik <- fit$loglik
        } else {
            fit <- fpot(x = y, threshold = 0, model = "gpd")
            est.y <- fit$estimate
            cov0.y <-  fit$var.cov
            colnames(cov0.y) <- rownames(cov0.y) <- names(est.y)
            logLik <- -fit$deviance / 2
        } 
        
    } else if(distname.y == "GPD") {
        
        ##=====================================================================
        ## Specific maximum likelihood for known shape (concave loglik),
        ## ML with 'evd::fpot' for the general case
        ##=====================================================================
        
        if (fixed.y["scale"]) {
            stop("parameter 'scale' can not be fixed alone in \"GPD\"")
        } else if (fixed.y["shape"]) {
            shape.hat <- fixed.par.y["shape"]
            fit <- fgpd1(x = y, shape = fixed.par.y["shape"])
            scale.hat <- fit$estimate
            est.y <- c(scale.hat, shape.hat)
            names(est.y) <- c("scale", "shape")
            cov0.y <- matrix(c(fit$cov, 0, 0, 0), nrow = 2L, ncol = 2L)
            colnames(cov0.y) <- rownames(cov0.y) <- names(est.y)
            logLik <- fit$loglik
        } else {
            fit <- fGPD(x = y)
            est.y <- fit$estimate
            cov0.y <-  fit$cov
            logLik <- fit$loglik
        } 
        
    } else if (distname.y == "gamma"){
        
        ##====================================================================
        ## Explicit maximum likelihood for known shape,
        ## ML with concentration with 'fgamma' for the general case
        ##====================================================================
        
        if (fixed.y["shape"]) {
            
            shape.hat <- fixed.par.y["shape"]
            if (shape.hat <= 0) stop("'fixed.par.y' must specify a 'shape' > 0")
            
            if (fixed.y["scale"]) {
                scale.hat <- fixed.par.y["scale"]
                cov0.y <- matrix(0, nrow = 2L, ncol = 2L)
                ymod <- y / scale.hat
                logLik <- n * ( -log(gamma(shape.hat)) - log(scale.hat) +
                               (shape.hat - 1) * mean(log(y)) - mean(ymod))
                
            } else {
                scale.hat <- mean(y)/shape.hat
                cov0.y <- matrix(0, nrow = 2L, ncol = 2L)
                cov0.y[2, 2] <- scale.hat * scale.hat / shape.hat / n
                logLik <- n * ( -log(gamma(shape.hat)) - log(scale.hat) +
                               (shape.hat - 1) * mean(log(y)) - shape.hat)
            }
            
            est.y <- c(shape = shape.hat, scale = scale.hat)
            names(est.y) <- parnames.y
            colnames(cov0.y) <- rownames(cov0.y) <- parnames.y
            
            
        } else if (fixed.y["scale"]) {
            stop("parameter 'scale' can not be fixed alone in gamma")
        } else {
            fit <- fgamma(x = y)
            est.y <- fit$estimate
            cov0.y <- fit$cov
            logLik <- fit$loglik
        }
        
        
    } else if (distname.y == "lomax"){
        
        ##======================================================================
        ## LOMAX estimation
        ## o Explicit maximum likelihood for the "known scale" case,
        ## o ML estimation (one parameter) for the "known shape" case.
        ## o ML with concentration with 'flomax' for the general case.
        ##======================================================================
        
        if (fixed.y["scale"]) {
            scale.hat <- fixed.par.y["scale"]
            ## we know that z follows an exponential distribution
            ## with rate = the shape parameter of the Lomax
            z <- log(1 + y / scale.hat)
            shape.hat <- 1 / mean(z)
            cov0.y <- matrix(0, nrow = 2L, ncol = 2L)
            cov0.y[1L, 1L] <- shape.hat*shape.hat / n
            logLik <- n * (log(shape.hat) - log(scale.hat) -
                     1.0 - 1.0 / shape.hat )
        } else if (fixed.y["shape"]) {
            shape.hat <- fixed.par.y["shape"]
            fit <- flomax1(x = y, shape = shape.hat)
            scale.hat <- fit$estimate
            cov0.y <- matrix(0, nrow = 2L, ncol = 2L)
            cov0.y[2L, 2L] <- fit$cov
            logLik <- fit$loglik
        } else {
            fit <- flomax(x = y, info.observed = FALSE)
            est <- fit$estimate
            shape.hat <- est["shape"]
            scale.hat <- est["scale"]
            cov0.y <- fit$cov
            logLik <- fit$loglik
        }
        ## CAUTION: here respect the order for Lomax!
        est.y <- c(shape.hat, scale.hat)
        names(est.y) <- parnames.y
        colnames(cov0.y) <- rownames(cov0.y) <- parnames.y
        
    } else if (distname.y == "maxlo"){
        
        if (fixed.y["scale"]) {
            scale.hat <- fixed.par.y["scale"]
            ## we know that z follows an exponential distribution
            ## with rate = the shape parameter of the Lomax
            if (any(y > scale.hat)) {
                stop("values of a maxlo distibution  must be < 'shape'")
            }
            z <- -log(1 - y/scale.hat)
            shape.hat <- 1 / mean(z)
            cov0.y <- matrix(0, nrow = 2L, ncol = 2L)
            cov0.y[1L, 1L] <- shape.hat*shape.hat / n
            logLik <- n * ( log(shape.hat) - log(scale.hat) - 1.0 + 1.0 / shape.hat )
        } else if (fixed.y["shape"]) {
            shape.hat <- fixed.par.y["shape"]
            fit <- fmaxlo1(x = y, shape = shape.hat)
            scale.hat <- fit$estimate
            cov0.y <- matrix(0, nrow = 2L, ncol = 2L)
            cov0.y[2L, 2L] <- fit$cov
            logLik <- fit$loglik
        } else {
            fit <- fmaxlo(x = y, info.observed = FALSE)
            est <- fit$estimate
            shape.hat <- est["shape"]
            scale.hat <- est["scale"]
            cov0.y <- fit$cov
            logLik <- fit$loglik
        }
        ## CAUTION: here respect the order for maxlo!
        est.y <- c(shape.hat, scale.hat)
        names(est.y) <- parnames.y
        colnames(cov0.y) <- rownames(cov0.y) <- parnames.y    
        
    }
    
    list(estimate = est.y,
         cov = cov0.y,
         logLik = logLik)
    
}
