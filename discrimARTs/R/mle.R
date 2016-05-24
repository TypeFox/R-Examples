mix.mle <- function(input, method='normal', 
    mix.prob=NULL, dist1.par1=NULL, dist1.par2=NULL, dist2.par1=NULL, dist2.par2=NULL, 
    lower=NULL, upper=NULL, 
    distlist=NULL, 
    optim.lower=rep(0,5), optim.upper=c(1, rep(NA,4)), ...) {
    ## use optim to find MLE value of 5 parameters 
    ## (mix.prob + (2parames *2 distributions))
    ## given the input measured trait vector and a method.
    ## mix.prob will be estimated from median(input) if not provided
    ## For method = normal, the remaining 4 distribution parameters will be estimated
    ## if none are provided 
    ## optionally, a list of 2 distribution pdf functions can be provided.
    ## Any additional options are passed to optim
    ## Return optimization object + original data and distlist

    ## pack distribution parameters into list for optim
    ## Order is important for bounds...
    ## Names are also important
    pars.init <- list(dist1.par1=dist1.par1, dist1.par2=dist1.par2, 
        dist2.par1=dist2.par1, dist2.par2=dist2.par2)
    ## Which parameters have *not* been supplied
    null.pars <-  sapply(pars.init, is.null)

    ## We can estimate all or none of the 4 distribution parameters
    if ( any(null.pars) && !all(null.pars) ) {
        stop('Distribution parameters  must either all be supplied, or none supplied (excluding mix.prob).')
    }

    ## We can only estimate for normal
    if ( all(null.pars) && method !='normal') {
        stop('Can only estimate distribution parameters for method="normal"')
    }
        
    if ( all(null.pars) && method=='normal') {
        message('Estimating initial conditions, assuming mix of normals')
        ## assume normal, estimate 
        ##??sd/2
        .sd <- sd(input)/sqrt(2)
        ## use as means of mixtures
        .quants <- quantile(input, probs=c(0.25, 0.75))
        names(.quants) <- NULL
        pars.init <- list(
            ## mean and sd of lower normal
            dist1.par1=.quants[1], dist1.par2=.sd,
            dist2.par1=.quants[2], dist2.par2=.sd
        )
    }

    ## Now add mix.prob to list, estimating if needed
    ## Needs to go first
    if (is.null(mix.prob)) {
        ## Estimate mixture from normal if absent
        ## add to pars.init later
        pars.init <- c(mix.prob=median(input), pars.init)
    } else {
        pars.init <- c(mix.prob=mix.prob, pars.init)
    }

    ## allow a user-supplied distlist
    if (!is.null(distlist)) {
        message('Ignoring method, using distlist')
        method <- 'user'
    } else {
        ## Use method, matching to internal distlist function via switch
        ## Match method to distlist function, call
        distlist <- switch( method, 
            normal = mix.distlist.norm(),
            facing.gamma = { 
                if (any(is.null(c(lower, upper)))) {
                    stop('lower and upper must be provided for method="facing.gamma"')
                }
                mix.distlist.facing.gamma( lower, upper)
            }
        )
        if (is.null(distlist)) { 
            ## No match from switch
            stop( sprintf( "Method %s not implemented", method))
        }
    }

    ## given an input vector of data, a vector of initial conditions, 
    ## and a list containing 2 PDF functions
    ## optimize mix.loglik over pars.init
    ## all other args are passed to optim
    ret <- optim(pars.init, mix.loglik, 
        .input=input,
        .distlist=distlist,
        hessian=T,
        ## 
        ## define sensible bounds
        ## 
        method='L-BFGS-B',
        lower=optim.lower,
        upper=optim.upper ,
        ...
    )
    ## Bind all 
    ret$method <- method
    ret$distlist <- distlist
    ret$pars.init <- pars.init
    ret$input <- input
    ## optim returns optimized values as par, rename for clarity
    ret$MLE.est <- ret$par
    ## optim returns negative logLik as value, rename for clarity 
    ret$neglogLik <- ret$value
    ## clean up
    ret$par <- NULL
    ret$value <- NULL
    ## Add bounds if available
    if (method=='facing.gamma') {
        ret$lower <- lower
        ret$upper <- upper
    }
    class(ret) <- 'discrimARTs'
    return(ret)
}
