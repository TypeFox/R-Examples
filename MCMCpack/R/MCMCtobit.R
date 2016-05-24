##########################################################################
## tobit regression model
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################



"MCMCtobit" <-
function(formula, data=NULL, below = 0.0, above = Inf,
           burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, ...) {

    # checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    if (!is.numeric(below) | !is.numeric(above)) {
        cat("Error: Censoring points must be numeric, which includes +-Inf.\n")
        stop("Please respecify and call ", calling.function(), " again.",
        call.=FALSE)    
    }
    if (below >= above) {
        cat("Error: Truncation points are logically inconsistent.\n")
        stop("Please respecify and call ", calling.function(), " again.",
        call.=FALSE)
    }
    
    # convert infinite values to finite approximations
    if(is.infinite(below)) below <- .Machine$double.xmax*-1
    if(is.infinite(above)) above <- .Machine$double.xmax
    
    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    # form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]    
    K <- ncol(X)  # number of covariates

    # starting values and priors
    beta.start <- coef.start(beta.start, K, formula, family=gaussian, data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    check.ig.prior(c0, d0)

    # define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, K+1)
    posterior <- NULL
    
    # call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCtobit", 
                     sample.nonconst=sample, Y=Y, X=X, below=as.double(below),
                     above=as.double(above),
                     burnin=as.integer(burnin), mcmc=as.integer(mcmc),
                     thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, c0=as.double(c0), d0=as.double(d0))
     
    # pull together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior,
                               names=c(xnames, "sigma2"),
                               title="MCMCtobit Posterior Sample")
    return(output)
  }
