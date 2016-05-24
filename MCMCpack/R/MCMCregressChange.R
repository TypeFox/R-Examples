#########################################################
## sample from the posterior distribution
## of a linear Gaussian model with multiple changepoints
## using linked C++ code in Scythe
##
## JHP 07/01/2007
## JHP 03/03/2009
#########################################################

"MCMCregressChange"<-
  function(formula, data=parent.frame(), m = 1, 
           b0 = 0, B0 = 0, c0 = 0.001, d0 = 0.001, sigma.mu = NA, sigma.var = NA, 
           a = NULL, b = NULL,
           mcmc = 1000, burnin = 1000,  thin = 1, verbose = 0, 
           seed = NA, beta.start = NA, P.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...){
    
    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    k <- ncol(X)    # number of covariates
    n <- length(y)
    ns <- m + 1 # number of states

    ## check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    cl <- match.call()
    nstore <- mcmc/thin    
    
    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    if(!is.na(seed)) set.seed(seed)

    ## prior
    mvn.prior <- form.mvn.prior(b0, B0, k)
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
       
    ## following MCMCregress, set chib as binary
    logmarglike <- loglik <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    if (m == 0){
      output <- MCMCregress(formula, burnin = burnin, mcmc = mcmc,
                            thin = thin, verbose = verbose,
                            b0 = b0, B0 = B0, c0 =c0, d0=d0,
                            marginal.likelihood = "Chib95")
      attr(output, "y") <- y
      attr(output, "m") <- m
    }
    else{
      if(k == 1){
        output <- MCMCresidualBreakAnalysis(y, data=data,  m = m, 
                                            b0 = b0, B0 = B0, c0 = c0, d0 = d0,
                                            a = a, b = b,
                                            burnin = burnin, mcmc = mcmc, thin = thin, verbose = verbose,
                                            seed = seed, beta.start = beta.start, P.start = P.start,
                                            marginal.likelihood = marginal.likelihood)
      }
      else{
        ## initial values
        Pstart  <-  check.P(P.start, m, a=a, b=b)
        A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)        
        betastart  <- beta.change.start(beta.start, ns, k, formula, family=gaussian, data)
        ols <- lm(y~X-1)
        Sigmastart <- rep(summary(ols)$sigma^2, ns)
        statestart <- sort(sample(1:ns, n, replace=T))
        
        ## call C++ code to draw sample
        posterior <- .C("MCMCregressChange", 
                        betaout = as.double(rep(0.0, nstore*ns*k)), 
                        Sigmaout = as.double(rep(0.0, nstore*ns)), 
                        Pout = as.double(rep(0.0, nstore*ns*ns)), 
                        psout = as.double(rep(0.0, n*ns)),
                        sout = as.double(rep(0.0, nstore*n)),
                        
                        Ydata = as.double(y),
                        Yrow = as.integer(nrow(y)),
                        Ycol = as.integer(ncol(y)),
                        Xdata = as.double(X),
                        Xrow = as.integer(nrow(X)),
                        Xcol = as.integer(ncol(X)),
                        
                        m = as.integer(m),
                        burnin = as.integer(burnin),           
                        mcmc = as.integer(mcmc), 
                        thin = as.integer(thin),
                        verbose = as.integer(verbose),
                        
                        lecuyer=as.integer(lecuyer), 
                        seedarray=as.integer(seed.array),
                        lecuyerstream=as.integer(lecuyer.stream),
                        
                        betastart = as.double(betastart),  
                        Sigmastart = as.double(Sigmastart),  
                        Pstart = as.double(Pstart),
                        statestart = as.integer(statestart),    
                        
                        a = as.double(a),
                        b = as.double(b),
                        b0data = as.double(b0),
                        B0data = as.double(B0), 
                        c0 = as.double(c0),
                        d0 = as.double(d0),
                        A0data = as.double(A0), 
                        logmarglikeholder = as.double(0.0),
                        loglikeholder = as.double(0.0),
                        chib = as.integer(chib))                
        
        ## get marginal likelihood if Chib95
        if (marginal.likelihood == "Chib95"){
          logmarglike <- posterior$logmarglikeholder
          loglike <- posterior$loglikeholder
        }
        
        ## pull together matrix and build MCMC object to return
        beta.holder <- matrix(posterior$betaout, nstore, ns*k)
        Sigma.holder <- matrix(posterior$Sigmaout, nstore, ns)
        P.holder    <- matrix(posterior$Pout, nstore, )
        s.holder    <- matrix(posterior$sout, nstore, )
        ps.holder   <- matrix(posterior$psout, n, )
        
        output1 <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
        varnames(output1)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c(xnames), "_regime", i, sep = "")
                                     })
        output2 <- mcmc(data=Sigma.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
        varnames(output2)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c("sigma2"), "_regime", i, sep = "")
                                     })
        output <- as.mcmc(cbind(output1, output2))
        
        attr(output, "title") <- "MCMCregressChange Posterior Sample"
        attr(output, "formula") <- formula
        attr(output, "y")       <- y
        attr(output, "X")       <- X
        attr(output, "m")       <- m
        attr(output, "call")    <- cl
        attr(output, "logmarglike") <- logmarglike
        attr(output, "loglik") <- loglik
        attr(output, "prob.state") <- ps.holder/nstore
        attr(output, "s.store") <- s.holder
      }
    }
    return(output)
      
  }## end of MCMC function
    
    
