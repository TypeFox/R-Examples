## sample from the posterior distribution
## of a binary model with multiple changepoints
## using linked C++ code in Scythe
##
## Written by JHP 07/01/2007
## Revised by JHP 07/16/2009

"MCMCbinaryChange"<-
  function(data,  m = 1, c0 = 1,  d0 = 1,  a = NULL, b = NULL, 
           burnin = 10000, mcmc = 10000, thin = 1, verbose = 0, 
           seed = NA, phi.start = NA, P.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...) {
    
    ## check iteration parameters
    check.mcmc.parameters(burnin, mcmc, thin)
    totiter <- mcmc + burnin
    cl <- match.call()
    
    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    if(!is.na(seed)) set.seed(seed)

    ## sample size
    y <- as.matrix(data)
    n <- nrow(y)
    ns <- m+1

    
    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    
    ## following MCMCregress, set chib as binary
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    nstore <- mcmc/thin    
    if (m == 0){
      b0 <- c0/(c0 + d0)
      B0 <- c0*d0/(c0 + d0)^2*(c0 + d0 + 1)
      output <- MCMCprobit(y~1, burnin = burnin, mcmc = mcmc,
                           thin = thin, verbose = verbose, b0 = b0, B0 = B0, 
                           marginal.likelihood = marginal.likelihood)
      attr(output, "y") <- y
    }
    else {
      ## prior for transition matrix
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)
      Pstart <- check.P(P.start, m=m, a=a, b=b)
      phistart <- check.theta(phi.start, ns, y, min=0, max=1)
    
    
      ## call C++ code to draw sample
      posterior <- .C("MCMCbinaryChange", 
                      phiout = as.double(rep(0.0, nstore*ns)), 
                      Pout = as.double(rep(0.0, nstore*ns*ns)), 
                      psout = as.double(rep(0.0, n*ns)),
                      sout = as.double(rep(0.0, nstore*n)), 
                      Ydata = as.double(y),
                      Yrow = as.integer(nrow(y)),
                      Ycol = as.integer(ncol(y)),
                      m = as.integer(m), 
                      burnin = as.integer(burnin),           
                      mcmc = as.integer(mcmc), 
                      thin = as.integer(thin),
                      verbose = as.integer(verbose),
                      lecuyer=as.integer(lecuyer), 
                      seedarray=as.integer(seed.array),
                      lecuyerstream=as.integer(lecuyer.stream),
                      phistart = as.double(phistart),  
                      Pstart = as.double(Pstart),    
                      a = as.double(a),
                      b = as.double(b),
                      c0 = as.double(c0),
                      d0 = as.double(d0), 
                      A0data = as.double(A0), 
                      logmarglikeholder = as.double(0.0),
                      chib = as.integer(chib))       
      
      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
      }
      
      ## pull together matrix and build MCMC object to return
      phi.holder <- matrix(posterior$phiout, nstore, )
      P.holder    <- matrix(posterior$Pout,  nstore, )
      s.holder    <- matrix(posterior$sout,  nstore, )
      ps.holder   <- matrix(posterior$psout, n, )
      
      output <- mcmc(data=phi.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output)  <- paste("phi.", 1:ns, sep = "")
      attr(output,"title") <- "MCMCbinaryChange Posterior Sample"
      attr(output, "y")       <- y
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      attr(output, "logmarglike") <- logmarglike
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "s.store") <- s.holder
    }
    return(output)
  }## end of MCMC function

