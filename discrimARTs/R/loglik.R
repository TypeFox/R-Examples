mix.loglik <- function(.pars, .input, .distlist, do.sum=TRUE) {
    ## optim expects pars vector is first argument
    ## Likelihood function for minimization by nlminb
    ## pars is the parameter vector  to be optimized 
    ## and (coerced to doubles) by nlminb
    ##
    ## optim passes along names of pars -- turn into named list
    ## and call parameters by name rather than position
    ## speed cost??
    parlist <- as.list(.pars)
    ## 
    ## likelihood of mixture 
    ## of given distributions for data
    #browser()
    lik <- with(parlist, 
            (1-mix.prob) * .distlist[[1]](.input, dist1.par1, dist1.par2) + 
            mix.prob     * .distlist[[2]](.input, dist2.par1, dist2.par2)
    )
    #browser()
    ## negative log likelihood
    ## (-2*loglik) is chi^2 distrib??
    ## ?? squared
    ret <- (-2*log(lik))
    ## return the total likelihood 
    ## or individual observation likelihood vector?
    if (do.sum) { 
        ret <- sum(ret)
    }
    return( ret )
}
