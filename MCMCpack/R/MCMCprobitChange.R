#########################################################
## 
## sample from the posterior distribution
## of a probit regression model with multiple changepoints
## 
## JHP 07/01/2007
## JHP 03/03/2009
#########################################################

"MCMCprobitChange"<-
  function(formula, data=parent.frame(),  m = 1, 
           burnin = 10000, mcmc = 10000, thin = 1, verbose = 0, 
           seed = NA, beta.start = NA, P.start = NA, 
           b0 = NULL, B0 = NULL, a = NULL, b = NULL,
           marginal.likelihood = c("none", "Chib95"), 
           ...){
    
    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    k <- ncol(X)   
    n <- length(y)
    ns <- m + 1

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

    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    if (m == 0){
      if (marginal.likelihood == "Chib95"){
        output <- MCMCprobit(formula=Y~X-1, burnin = burnin, mcmc = mcmc,
                             thin = thin, verbose = verbose,
                             b0 = b0, B0 = B0,
                             marginal.likelihood = "Laplace")
        cat("\n Chib95 method is not yet available for m = 0. Laplace method is used instead.")
      }
      else {
        output <- MCMCprobit(formula=Y~X-1, burnin = burnin, mcmc = mcmc,
                             thin = thin, verbose = verbose,
                             b0 = b0, B0 = B0)
      }
      attr(output, "y") <- y
    }
    else{
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)  
      Pstart  <-  check.P(P.start, m, a=a, b=b)
      betastart  <- beta.change.start(beta.start, ns, k, formula, family=binomial, data)
      ## call C++ code to draw sample
      posterior <- .C("MCMCprobitChange", 
                      betaout = as.double(rep(0.0, nstore*ns*k)), 
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
                      Pstart = as.double(Pstart),    
                      
                      a = as.double(a),
                      b = as.double(b),
                      b0data = as.double(b0),
                      B0data = as.double(B0), 
                      A0data = as.double(A0), 
                      logmarglikeholder = as.double(0.0),
                      loglikeholder = as.double(0.0),
                      chib = as.integer(chib))                
      
      ## get marginal likelihood if Chib95
      if (chib==1){
        logmarglike <- posterior$logmarglikeholder
        loglike <- posterior$loglikeholder
      }
      else{
        logmarglike <- loglike <- 0
      }
      
      ## pull together matrix and build MCMC object to return
      beta.holder <- matrix(posterior$betaout, nstore, ns*k)
      P.holder    <- matrix(posterior$Pout, nstore, )
      s.holder    <- matrix(posterior$sout, nstore, )
      ps.holder   <- matrix(posterior$psout, n, )
      
      output <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output)  <- sapply(c(1:ns),
                                  function(i){
                                    paste(c(xnames), "_regime", i, sep = "")
                                  })
      attr(output, "title") <- "MCMCprobitChange Posterior Sample"
      attr(output, "formula") <- formula
      attr(output, "y")       <- y
      attr(output, "X")       <- X
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglike") <- loglike
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "s.store") <- s.holder
    }
    return(output)

}## end of MCMC function

