rvalueBoot <- function(object, statistic = median, R, type = "nonparametric") {
    if(object$aux$family == "gaussian") {
        estimate <- object$aux$unsorted$MLE
        nuisance <- object$aux$unsorted$SE
    }
    else if(object$aux$family == "binomial")  {
        estimate <- object$aux$unsorted$xx
        nuisance <- object$aux$unsorted$nn
    }
    else if(object$aux$family == "poisson") {
        estimate <- object$aux$unsorted$xx
        nuisance <- object$aux$unsorted$eta
    }
    else if(object$aux$family == "tdist") {
        estimate <- object$aux$unsorted$MLE
        nuisance <- object$aux$unsorted$SE
    }
    nunits <- length(estimate)

    RvalStore <- matrix(nrow = nunits, ncol = R+1)
    RvalStore[,1] <- object$rvalues
    nn <- seq_len(nunits)
    dd <- ddtmp <- cbind(estimate,nuisance)
    if(type=="nonparametric") {
        if(object$aux$prior=="conjugate") {
            for(i in 1:R) {
               ### need to add parametric bootstrap.
               subsamp.index <- sample(nn, size = nunits, replace = TRUE)
               ddtmp[,1] <- estimate[subsamp.index]
               ddtmp[,2] <- nuisance[subsamp.index]
               new.hypers <- HyperEstimate(ddtmp[,1],ddtmp[,2], family = object$aux$family)
               temp <- rvalues(dd,family = object$aux$family, prior = object$aux$prior,
                              hypers = new.hypers, alpha.grid = object$aux$alpha.grid,
                              smooth = object$aux$smooth)
               ### need to add other parameters to rvalues function (e.g. smooth)
               RvalStore[,i+1] <- temp$rvalues
            }
        }
        else if(object$aux$prior=="nonparametric") {
            for(i in 1:R) {
                subsamp.index <- sample(nn, size = nunits, replace = TRUE)
                ddtmp[,1] <- estimate[subsamp.index]
                ddtmp[,2] <- nuisance[subsamp.index]
                
                if(is.null(object$aux$df)) {
                    npfit <- npmle(ddtmp,family = object$aux$family)
                }
                else {
                    npfit <- npmle(ddtmp,family = tdist(df=object$aux$df))
                }
                lb <- min(npfit$support) - .01
                theta.alpha <- ThetaQuantiles(npfit$Fhat, alpha.grid = object$aux$alpha.grid, 
                                              lbd = lb, ubd = max(npfit$support))
             
                theta.probs <- -diff(c(1,npfit$Fhat(theta.alpha)))
                tmp <- NPagrid(estimate = dd[,1], nuisance = dd[,2], theta.alpha, 
                               theta.probs, alpha.grid = object$aux$alpha.grid, 
                               smooth = object$aux$smooth, family = object$aux$family,
                               df = object$aux$df)
                RvalStore[,i+1] <- tmp$rvalues
            }
        }
    }
    else if(type=="parametric") {
        ### The parametric bootstrap refers to generating replicate data sets in the following way:
        ### Generate \theta_i ~ pi(\theta| hyper_est). Then, generate X_i ~ p(x|\theta_i,\eta_i)
         
        if(object$aux$prior=="conjugate") {
            for(i in 1:R){
                switch(object$aux$family,
                gaussian={
                   theta <- rnorm(nunits,mean=object$aux$hypers[1],sd=sqrt(object$aux$hypers[2]))
                   new.estimate <- theta + nuisance*rnorm(nunits)
                },
                poisson={
                   theta <- rgamma(nunits,shape=object$aux$hypers[1],rate=object$aux$hypers[2])
                   new.estimate <- rpois(nunits,lambda=theta*nuisance)
                },
                binomial={
                   theta <- rbeta(nunits,shape1=object$aux$hypers[1],shape2=object$aux$hypers[2])
                   new.estimate <- rbinom(nunits,prob=theta,size=nuisance)
                },
                )
                
                new.hypers <- HyperEstimate(new.estimate,nuisance,family=object$aux$family)
                tmp <- rvalues(dd,family = object$aux$family, prior = object$aux$prior,
                               hypers = new.hypers,alpha.grid = object$aux$alpha.grid,
                               smooth = object$aux$smooth)
                RvalStore[,i+1] <- tmp$rvalues
            }
        }
        else if(object$aux$prior=="nonparametric"){
            for(i in 1:R){
                ### suppress the alias method warning (used when there are 
                ### more than 250 reasonably probable values)
                theta <- suppressWarnings(sample(object$aux$support, size = nunits, replace = TRUE, prob = object$aux$mix.prop))
                switch(object$aux$family,
                gaussian={
                   ddtmp[,1] <- theta + nuisance*rnorm(nunits)
                   npfit <- npmle(ddtmp,family = object$aux$family)
                },
                poisson={
                   ddtmp[,1] <- rpois(nunits,lambda=theta*nuisance)
                   npfit <- npmle(ddtmp,family = poisson)
                },
                binomial={
                   ddtmp[,1] <- rbinom(nunits,prob=theta,size=nuisance)
                   npfit <- npmle(ddtmp,family = object$aux$family)
                },
                tdist={
                   ddtmp[,1] <- theta + nuisance*rt(nunits,df=object$aux$df)
                   npfit <- npmle(ddtmp,family = tdist(df=object$aux$df))
                },
                )
                lb <- min(npfit$support) - .01
                theta.alpha <- ThetaQuantiles(npfit$Fhat, alpha.grid = object$aux$alpha.grid, 
                                              lbd = lb, ubd = max(npfit$support))
             
                theta.probs <- -diff(c(1,npfit$Fhat(theta.alpha)))
                tmp <- NPagrid(estimate = dd[,1], nuisance = dd[,2], theta.alpha, 
                               theta.probs, alpha.grid = object$aux$alpha.grid, 
                               smooth = object$aux$smooth, family = object$aux$family,
                               df=object$aux$df)
                
                RvalStore[,i+1] <- tmp$rvalues
            }
        }
    }
    else{ 
        stop("Bootstrap type must be either nonparametric or parametric.")
    }
    ans <- list()
    ans$rval.repmat <- RvalStore
    ans$rval.boot <- apply(RvalStore,1,FUN=statistic)
    return(ans)
}
### Testing commit on 12/16/15
#### To do: add control parameters to the npmle 
