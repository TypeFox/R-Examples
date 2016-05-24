#########################################################
## residual break analysis
## JHP 07/01/2007
## JHP 03/03/2009
#########################################################

"MCMCresidualBreakAnalysis"<-
  function(resid, m = 1, 
           b0 = 0, B0 = 0.001, c0 = 0.1, d0 = 0.1, a = NULL, b = NULL,
           mcmc = 1000, burnin = 1000,  thin = 1, verbose = 0, 
           seed = NA, beta.start = NA, P.start = NA,
           marginal.likelihood = c("none", "Chib95"), ...){
    
    ## form response and model matrices
    y <- as.matrix(resid, ,1)
    n <- nrow(y)
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

   ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    
    ## following MCMCregress, set chib as binary
    logmarglike <- loglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }
    
    ## initial values
    if(m == 0){
      output <- MCMCregress(y~1, mcmc=mcmc, burnin=burnin, verbose=verbose, thin=thin,
                            b0=b0, B0=B0, c0=c0, d0=d0, seed = seed,
                            beta.start=beta.start, marginal.likelihood = marginal.likelihood)
      attr(output, "y") <- y
    }
    else{
      ## prior
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)  
      Pstart  <-  check.P(P.start, m, a=a, b=b)
      betastart  <- rep(mean(y), ns)
      Sigmastart <- rep(var(y), ns)
      statestart <- sort(sample(1:ns, n, replace=T))
      
      ## call C++ code to draw sample
      posterior <- .C("MCMCresidualBreakAnalysis", 
                      betaout = as.double(rep(0.0, nstore*ns)), 
                      Sigmaout = as.double(rep(0.0, nstore*ns)), 
                      psout = as.double(rep(0.0, n*ns)),
                      
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
      beta.holder <- matrix(posterior$betaout, nstore, ns)
      Sigma.holder <- matrix(posterior$Sigmaout, nstore, ns)
      ps.holder   <- matrix(posterior$psout, n, )
      
      output1 <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output1)  <- sapply(c(1:ns),
                                   function(i){
                                     paste("mu_regime", i, sep = "")
                                   })
      output2 <- mcmc(data=Sigma.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output2)  <- sapply(c(1:ns),
                                   function(i){
                                   paste("sigma2_regime", i, sep = "")
                                 })
      output <- as.mcmc(cbind(output1, output2))
      
      attr(output, "title") <- "MCMCresidualBreakAnalysis Posterior Sample"
      attr(output, "formula") <- formula
      attr(output, "y")       <- y
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglike") <- loglike
    }
    return(output)
    
}## end of MCMC function

