##########################################################################
## sample from the posterior distribution of a probit
## model in R using linked C++ code in Scythe
##
## ADM and KQ 5/21/2002
## Modified to meet new developer specification 7/26/2004 KQ
## Modified for new Scythe and rngs 7/26/2004 KQ
## Modified to handle marginal likelihood calculation 1/27/2006 KQ
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"MCMCprobit" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin = 1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, bayes.resid=FALSE,
           marginal.likelihood = c("none", "Laplace", "Chib95"), ...) {
    
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
    beta.start <- coef.start(beta.start, K, formula,
                             family=binomial(link="probit"), data)
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

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

    
    ## residuals setup
    resvec <- NULL
    if (is.logical(bayes.resid) && bayes.resid==TRUE){
      resvec <- matrix(1:length(Y), length(Y), 1)
    }
    else if (!is.logical(bayes.resid)){
      resvec <- matrix(bayes.resid, length(bayes.resid), 1)
      if (min(resvec %in% 1:length(Y)) == 0){
        cat("Elements of bayes.resid are not valid row numbers.\n")
        stop("Check data and call MCMCprobit() again.\n") 
      }
    }
    
    ## y \in {0, 1} error checking
    if (sum(Y!=0 & Y!=1) > 0) {
      cat("Elements of Y equal to something other than 0 or 1.\n")
      stop("Check data and call MCMCprobit() again.\n") 
    }

    ## if Chib95 is true
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    posterior <- NULL

    if (is.null(resvec)){
      ## define holder for posterior density sample
      sample <- matrix(data=0, mcmc/thin, dim(X)[2] )
  
      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCprobit",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose), betastart=beta.start,
                       b0=b0, B0=B0, logmarglikeholder.nonconst = as.double(0.0),
                       chib = as.integer(chib)) 
      
      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
      }
      ## marginal likelihood calculation if Laplace
      if (marginal.likelihood == "Laplace"){
        theta.start <- beta.start
        optim.out <- optim(theta.start, logpost.probit, method="BFGS",
                           control=list(fnscale=-1),
                           hessian=TRUE, y=Y, X=X, b0=b0, B0=B0)
        
        theta.tilde <- optim.out$par
        beta.tilde <- theta.tilde[1:K]
        
        Sigma.tilde <- solve(-1*optim.out$hessian)
        
        logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
          log(sqrt(det(Sigma.tilde))) + 
            logpost.probit(theta.tilde, Y, X, b0, B0)
        
      }
      
      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCprobit Posterior Sample",
                                 y=Y, call=cl,
                                 logmarglike=logmarglike)

    }
    else{
      # define holder for posterior density sample
      sample <- matrix(data=0, mcmc/thin, dim(X)[2]+length(resvec) )
      
      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCprobitres",
                       sample.nonconst=sample, Y=Y, X=X, resvec=resvec,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose), betastart=beta.start,
                       b0=b0, B0=B0,  logmarglikeholder.nonconst= as.double(0.0),
                       chib = as.integer(chib)) 
      
      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
      }
      ## marginal likelihood calculation if Laplace
      if (marginal.likelihood == "Laplace"){
        theta.start <- beta.start
        optim.out <- optim(theta.start, logpost.probit, method="BFGS",
                           control=list(fnscale=-1),
                           hessian=TRUE, y=Y, X=X, b0=b0, B0=B0)
        
        theta.tilde <- optim.out$par
        beta.tilde <- theta.tilde[1:K]
        
        Sigma.tilde <- solve(-1*optim.out$hessian)
        
        logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
          log(sqrt(det(Sigma.tilde))) + 
            logpost.probit(theta.tilde, Y, X, b0, B0)
        
      }
      
      ## put together matrix and build MCMC object to return
      xnames <- c(xnames, paste("epsilonstar", as.character(resvec), sep="") )
      
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCprobit Posterior Sample",
                                 y=Y, call=cl, logmarglike=logmarglike)
      
    }
    return(output)
  
  }


