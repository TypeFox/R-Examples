################################
## Poisson Changepoint Model
##
## 07/14/2009 Jong Hee Park
################################

"MCMCpoissonChange"<-
  function(formula, data = parent.frame(), m = 1,
           b0 = 0, B0 = 1, a = NULL, b = NULL, c0 = NA, d0 = NA,
           lambda.mu = NA, lambda.var = NA, 
           burnin = 1000, mcmc = 1000, thin = 1, verbose = 0, 
           seed = NA, beta.start = NA, P.start = NA, ## offset = NA,   
           marginal.likelihood = c("none", "Chib95"), ...) {
    
    ## form response and model matrices
    holder <- parse.formula(formula, data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    k <- ncol(X)    
    n <- length(y)
    n.arrival<-  y + 1
    NT      <-  max(n.arrival)   
    tot.comp <-  n + sum(y)      
    ns      <- m + 1 
    
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
        
    if (k==1){
      if (!is.na(lambda.mu) && !is.na(lambda.var)) {
        c0 <- lambda.mu^2/lambda.var
        d0 <- lambda.mu/lambda.var
      }
      if ((is.na(c0)||is.na(d0))&&((is.na(lambda.mu)||is.na(lambda.var)))){
        stop("You have to provide a prior for lambda (c0 and d0 or lambda.mu and lambda.var) when there is no covariate.\n")
      }
    }
    else{
      c0 <- d0 <- 0
      mvn.prior <- form.mvn.prior(b0, B0, k)
      b0 <- mvn.prior[[1]]
      B0 <- mvn.prior[[2]]
    }
    
    
    
    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }
    if (m == 0){
      if (marginal.likelihood == "Chib95"){
        if (is.na(b0)||is.na(B0))
          stop("You have to have a prior for beta (b0 and B0) when m = 0.\n")
        else{
          output <- MCMCpoisson(formula, burnin = burnin, mcmc = mcmc,
                                thin = thin, verbose = verbose,
                                b0 = b0, B0 = B0,
                                marginal.likelihood = "Laplace")
          cat("\n Chib95 method is not yet available for m = 0. Laplace method is used instead.")
        }
      }
      else {
        output <- MCMCpoisson(formula, burnin = burnin, mcmc = mcmc,
                              thin = thin, verbose = verbose,
                              b0 = b0, B0 = B0)
      }
    }
    else {
      ## prior
      A0 <- trans.mat.prior(m=m, n=n, a=a, b=b)  

      ## get initial values of tau from observed y
      Pstart <- check.P(P.start, m, a=a, b=b)
      betastart  <- beta.change.start(beta.start, ns, k, formula, family=poisson, data)
      if (k == 1){
        betastart <- exp(betastart)
      }
      taustart <- tau.initial(y, tot.comp)
      componentstart  <-  round(runif(tot.comp, 1, 5))
      ## if(is.na(offset)){
      ##   logoffset <- rep(0, length(y))
      ## }
      ## else{
      ##    if(length(offset) == length(y)){
      ##      logoffset <- log(offset)
      ##   }
      ##   else{
      ##     stop("\n The length of offset is not matched with y.")
      ##   }
      ##  }
      ##  print(offset)
      
      ## normal mixture weights
      wr  <-  c(0.2294, 0.2590, 0.2480, 0.1525, 0.0472)
      mr  <-  c(0.0982, -1.5320, -0.7433, 0.8303, -3.1428)
      sr  <-  sqrt(c(0.2401, 1.1872, 0.3782, 0.1920, 3.2375))
      
      ## call C++ code to draw sample
      posterior <- .C("MCMCpoissonChange", 
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
                      ## logoffset = as.double(logoffset), 
                      m = as.integer(m), 
                      burnin = as.integer(burnin),           
                      mcmc = as.integer(mcmc), 
                      thin = as.integer(thin),
                      verbose = as.integer(verbose),                     
                      betastart = as.double(betastart),  
                      Pstart = as.double(Pstart),    
                      taustart = as.double(taustart),    
                      componentstart = as.double(componentstart),    
                      a = as.double(a),
                      b = as.double(b),
                      c0 = as.double(c0),
                      d0 = as.double(d0),
                      lecuyer = as.integer(lecuyer),
                      seedarray = as.integer(seed.array),
                      lecuyerstream = as.integer(lecuyer.stream),
                      b0data = as.double(b0),
                      B0data = as.double(B0), 
                      A0data = as.double(A0), 
                      logmarglikeholder = as.double(0.0),
                      loglikeholder = as.double(0.0),
                      wrin = as.double(wr),
                      mrin = as.double(mr),
                      srin = as.double(sr),
                      chib = as.integer(chib), 
                      
                      PACKAGE="MCMCpack")                
      
      ## get marginal likelihood if Chib95
      if (marginal.likelihood == "Chib95"){
        logmarglike <- posterior$logmarglikeholder
        loglike <- posterior$loglikeholder
      }
      else {
        logmarglike <- 0
        loglike <- 0
      }
 
      
      ## pull together matrix and build MCMC object to return
      beta.holder <- matrix(posterior$betaout, nstore, )
      P.holder    <- matrix(posterior$Pout, nstore, )
      s.holder    <- matrix(posterior$sout, nstore, )
      ps.holder   <- matrix(posterior$psout, n, )
      
      output <- mcmc(data=beta.holder, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output)  <- sapply(c(1:ns), function(i) { paste(xnames, "_regime", i, sep = "")})
      attr(output, "title") <- "MCMCpoissonChange Posterior Sample"
      attr(output, "formula") <- formula
      attr(output, "y")       <- y
      attr(output, "X")       <- X
      attr(output, "m")       <- m
      attr(output, "call")    <- cl
      attr(output, "logmarglike") <- logmarglike
      attr(output, "loglike") <- loglike
      attr(output, "prob.state") <- ps.holder/nstore
      attr(output, "s.store") <- s.holder
      attr(output, "P.store") <- P.holder
    }
    return(output)
      
  }
