##########################################################################
## sample from the posterior distribution of a logistic regression
## model in R using linked C++ code in Scythe
##
## KQ 1/23/2003
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Modified to meet new developer specification 7/15/2004 KQ
## Modified for new Scythe and rngs 7/25/2004 KQ
## note: B0 is now a precision
## Modified to allow user-specified prior density 8/17/2005 KQ
## Modified to handle marginal likelihood calculation 1/27/2006 KQ
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"MCMClogit" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin=1, tune=1.1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, user.prior.density=NULL, logfun=TRUE,
           marginal.likelihood = c("none", "Laplace"), ...) {

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
    beta.start <- coef.start(beta.start, K, formula, family=binomial, data)
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

    
    ## setup the environment so that fun can see the things passed as ...
    userfun <- function(ttt) user.prior.density(ttt, ...)
    my.env <- environment(fun = userfun)
    
    
    
    ## check to make sure user.prior.density returns a numeric scalar and
    ## starting values have positive prior mass
    if (is.function(user.prior.density)){
      funval <- userfun(beta.start)
      if (length(funval) != 1){
        cat("user.prior.density does not return a scalar.\n")
        stop("Respecify and call MCMClogit() again. \n")       
      }
      if (!is.numeric(funval)){
        cat("user.prior.density does not return a numeric value.\n")
        stop("Respecify and call MCMClogit() again. \n")       
      }

      if (identical(funval,  Inf)){
        cat("user.prior.density(beta.start) == Inf.\n")
        stop("Respecify and call MCMClogit() again. \n")       
      }
      
      if (logfun){
        if (identical(funval, -Inf)){
          cat("user.prior.density(beta.start) == -Inf.\n")
          stop("Respecify and call MCMClogit() again. \n")       
        }
      }
      else{
        if (funval <= 0){
          cat("user.prior.density(beta.start) <= 0.\n")
          stop("Respecify and call MCMClogit() again. \n")       
        }        
      }   
    }
    else if (!is.null(user.prior.density)){
      cat("user.prior.density is neither a NULL nor a function.\n")
      stop("Respecify and call MCMClogit() again. \n")       
    }
        

       
    ## form the tuning parameter
    tune <- vector.tune(tune, K)
    V <- vcov(glm(formula=formula, data=data, family=binomial))
      
    ## y \in {0, 1} error checking
    if (sum(Y!=0 & Y!=1) > 0) {
      cat("Elements of Y equal to something other than 0 or 1.\n")
      stop("Check data and call MCMClogit() again. \n") 
    }
   
   

    propvar <- tune %*% V %*% tune

    posterior <- NULL    
    ## call C++ code to draw sample
    if (is.null(user.prior.density)){
      ## define holder for posterior density sample
      sample <- matrix(data=0, mcmc/thin, dim(X)[2] )
      
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMClogit",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       tune=tune, lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose), betastart=beta.start,
                       b0=b0, B0=B0, V=V) 

      ## marginal likelihood calculation if Laplace
      if (marginal.likelihood == "Laplace"){
        theta.start <- beta.start
        optim.out <- optim(theta.start, logpost.logit, method="BFGS",
                           control=list(fnscale=-1),
                           hessian=TRUE, y=Y, X=X, b0=b0, B0=B0)
        
        theta.tilde <- optim.out$par
        beta.tilde <- theta.tilde[1:K]
        
        Sigma.tilde <- solve(-1*optim.out$hessian)
        
        logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
          log(sqrt(det(Sigma.tilde))) + 
            logpost.logit(theta.tilde, Y, X, b0, B0)
        
      }
      
      
      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMClogit Posterior Sample",
                                 y=Y, call=cl, logmarglike=logmarglike)
      
    }
    else {
      sample <- .Call("MCMClogituserprior_cc",
                      userfun, as.integer(Y), as.matrix(X),
                      as.double(beta.start),
                      my.env, as.integer(burnin), as.integer(mcmc),
                      as.integer(thin),
                      as.integer(verbose),
                      lecuyer=as.integer(lecuyer), 
                      seedarray=as.integer(seed.array),
                      lecuyerstream=as.integer(lecuyer.stream),
                      as.logical(logfun),
                      as.matrix(propvar),
                      PACKAGE="MCMCpack")

      
      ## marginal likelihood calculation if Laplace
      if (marginal.likelihood == "Laplace"){
        theta.start <- beta.start
        optim.out <- optim(theta.start, logpost.logit.userprior, method="BFGS",
                           control=list(fnscale=-1),
                           hessian=TRUE, y=Y, X=X, userfun=userfun,
                           logfun=logfun, my.env=my.env)
        
        theta.tilde <- optim.out$par
        beta.tilde <- theta.tilde[1:K]
        
        Sigma.tilde <- solve(-1*optim.out$hessian)
        
        logmarglike <- (length(theta.tilde)/2)*log(2*pi) +
          log(sqrt(det(Sigma.tilde))) + 
            logpost.logit.userprior(theta.tilde, Y, X, userfun=userfun,
                                    logfun=logfun, my.env=my.env)
        
      }

      
      output <- mcmc(data=sample, start=burnin+1,
                     end=burnin+mcmc, thin=thin)
      varnames(output) <- as.list(xnames)
      attr(output, "title") <- "MCMClogit Posterior Sample"
      attr(output, "y") <- Y
      attr(output, "call") <- cl
      attr(output, "logmarglike") <- logmarglike
    }
    
    return(output)    
  }

##########################################################################


