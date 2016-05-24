ST.CARlinear <- function(formula, family, data=NULL,  trials=NULL, W, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.mean.alpha=NULL, prior.var.alpha=NULL, prior.nu2=NULL, prior.tau2=NULL, fix.rho.slo=FALSE, rho.slo=NULL, fix.rho.int=FALSE, rho.int=NULL, verbose=TRUE)
{
    ## This is a wrapper function for the following three functions.
    ## binomial.CARanova
    ## gaussian.CARanova
    ## poisson.CARanova
    if(is.null(family)) stop("the family argument is missing", call.=FALSE)
    
    #### Run the appropriate model according to the family arugment
    if(family=="binomial")
    {
        if(is.null(trials)) stop("a binomial model was specified but the trials arugment was not specified", call.=FALSE)
        model <- binomial.CARlinear(formula=formula, data=data,  trials=trials, W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.alpha=prior.mean.alpha, prior.var.alpha=prior.var.alpha, prior.tau2=prior.tau2, fix.rho.slo=fix.rho.slo, rho.slo=rho.slo, fix.rho.int=fix.rho.int, rho.int=rho.int, verbose=verbose)
    }else if(family=="gaussian")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- gaussian.CARlinear(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.alpha=prior.mean.alpha, prior.var.alpha=prior.var.alpha, prior.nu2=prior.nu2, prior.tau2=prior.tau2, fix.rho.slo=fix.rho.slo, rho.slo=rho.slo, fix.rho.int=fix.rho.int, rho.int=rho.int, verbose=verbose)          
    }else if(family=="poisson")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- poisson.CARlinear(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.mean.alpha=prior.mean.alpha, prior.var.alpha=prior.var.alpha, prior.tau2=prior.tau2, fix.rho.slo=fix.rho.slo, rho.slo=rho.slo, fix.rho.int=fix.rho.int, rho.int=rho.int, verbose=verbose)          
    }else
    {
        stop("the family arugment is not one of `binomial', `gaussian' or `poisson'.", call.=FALSE)     
    }
    return(model)     
}