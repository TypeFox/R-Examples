ST.CARanova <- function(formula, family, data=NULL,  trials=NULL, W, interaction=TRUE, burnin, n.sample, thin=1, prior.mean.beta=NULL, prior.var.beta=NULL, prior.nu2=NULL, prior.tau2=NULL, fix.rho.S=FALSE, rho.S=NULL, fix.rho.T=FALSE, rho.T=NULL, verbose=TRUE)
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
        model <- binomial.CARanova(formula=formula, data=data,  trials=trials, W=W, interaction=interaction, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, fix.rho.S=fix.rho.S, rho.S=rho.S, fix.rho.T=fix.rho.T, rho.T=rho.T, verbose=verbose)
    }else if(family=="gaussian")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- gaussian.CARanova(formula=formula, data=data,  W=W, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.nu2=prior.nu2, prior.tau2=prior.tau2, fix.rho.S=fix.rho.S, rho.S=rho.S, fix.rho.T=fix.rho.T, rho.T=rho.T, verbose=verbose)          
    }else if(family=="poisson")
    {
        if(!is.null(trials)) stop("you do not need a trials arugment as a binomial model was not specified", call.=FALSE)
        model <- poisson.CARanova(formula=formula, data=data,  W=W, interaction=interaction, burnin=burnin, n.sample=n.sample, thin=thin, prior.mean.beta=prior.mean.beta, prior.var.beta=prior.var.beta, prior.tau2=prior.tau2, fix.rho.S=fix.rho.S, rho.S=rho.S, fix.rho.T=fix.rho.T, rho.T=rho.T, verbose=verbose)          
    }else
    {
        stop("the family arugment is not one of `binomial', `gaussian' or `poisson'.", call.=FALSE)     
    }
    return(model)     
}