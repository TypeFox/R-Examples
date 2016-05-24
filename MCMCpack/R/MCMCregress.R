##########################################################################
## MCMCregress.R samples from the posterior distribution of a Gaussian
## linear regression model in R using linked C++ code in Scythe
##
## Original written by ADM and KQ 5/21/2002
## Updated with helper functions ADM 5/28/2004
## Modified to meet new developer specification 6/18/2004 KQ
## Modified for new Scythe and rngs 7/22/2004 ADM
## Modified to handle marginal likelihood calculation 1/26/2006 KQ
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

"MCMCregress" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, sigma.mu = NA, sigma.var = NA, 
           marginal.likelihood = c("none", "Laplace", "Chib95"),
           ...) {
    
    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    
    cl <- match.call()
        
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]    
    K <- ncol(X)  # number of covariates
        
    ## starting values and priors
    beta.start <- coef.start(beta.start, K, formula, family=gaussian, data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    
    if (is.na(sigma.mu)|is.na(sigma.var)) {
      check.ig.prior(c0, d0)
    }
    else {
      c0 <- 4 + 2 *(sigma.mu^2/sigma.var)
      d0 <- 2*sigma.mu *(c0/2 - 1)
    }
      
    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    B0.eigenvalues <- eigen(B0)$values
    if (min(B0.eigenvalues) < 0){
      stop("B0 is not positive semi-definite.\nPlease respecify and call again.\n")
    }
    if (isTRUE(all.equal(min(B0.eigenvalues), 0))){
      if (marginal.likelihood != "none"){
        warning("Cannot calculate marginal likelihood with improper prior\n")
        marginal.likelihood <- "none"
      }
    }
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    
    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, K+1)
    posterior <- NULL 
    
    ## call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCregress", 
                     sample.nonconst=sample, Y=Y, X=X, 
                     burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer), 
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), betastart=beta.start,
                     b0=b0, B0=B0, c0=as.double(c0), d0=as.double(d0),
                     logmarglikeholder.nonconst=as.double(0.0),
                     chib=as.integer(chib))

    ## get marginal likelihood if Chib95
    if (marginal.likelihood == "Chib95"){
      logmarglike <- posterior$logmarglikeholder
    }
    ## marginal likelihood calculation if Laplace
    if (marginal.likelihood == "Laplace"){
      theta.start <- c(beta.start, log(0.5*var(Y)))
      optim.out <- optim(theta.start, logpost.regress, method="BFGS",
                         control=list(fnscale=-1),
                         hessian=TRUE, y=Y, X=X, b0=b0, B0=B0, c0=c0, d0=d0)

      theta.tilde <- optim.out$par
      beta.tilde <- theta.tilde[1:K]
      sigma2.tilde <- exp(theta.tilde[K+1])
      
      Sigma.tilde <- solve(-1*optim.out$hessian)
      
      logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
        log(sqrt(det(Sigma.tilde))) + 
          logpost.regress(theta.tilde, Y, X, b0, B0, c0, d0)
      
    }
    
    ## pull together matrix and build MCMC object to return
    output <- form.mcmc.object(posterior,
                               names=c(xnames, "sigma2"),
                               title="MCMCregress Posterior Sample",
                               y=Y, call=cl,
                               logmarglike=logmarglike
                               )
    return(output)
  }
