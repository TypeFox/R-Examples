#########################################################
## 
## sample from the posterior distribution
## of ordinal probit changepoint regression model
## using a linear Gaussian approximation
##
## JHP 07/01/2007
## JHP 03/03/2009
## JHP 09/08/2010
#########################################################

"MCMCoprobitChange"<-
  function(formula, data=parent.frame(),  m = 1, 
           burnin = 1000, mcmc = 1000, thin = 1, tune = NA, verbose = 0, 
           seed = NA, beta.start = NA, gamma.start = NA, P.start = NA, 
           b0 = NULL, B0 = NULL, a = NULL, b = NULL,
           marginal.likelihood = c("none", "Chib95"), 
           gamma.fixed=0, ...){

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    cl <- match.call()
    nstore <- mcmc/thin
      
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    totiter <- mcmc+burnin
    holder <- parse.formula(formula, data=data)
    y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    K <- ncol(X)
    Y <- factor(y, ordered = TRUE)
    ncat <- nlevels(Y) 
    cat <- levels(Y)
    ns <- m + 1
    N <- nrow(X)
    gk <- ncat + 1

    if(sum(is.na(tune))==1) {
      stop("Please specify a tune parameter and call MCMCoprobitChange() again.\n") 
    }
    else if (length(tune)==1){
      tune <- rep(tune, ns)
    }
    else if(length(tune)>1&length(tune)<ns){
      tune <- rep(tune[1], ns)
      cat("The first element of tune is repeated to make it conformable to the number of states.\n")
    }
    else{
     
    }
               
    xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if (xint > 0) {
      new.X <- X[, -xint, drop = FALSE]
    }
    else
      warning("An intercept is needed and assumed in MCMCoprobitChange()\n.")
    if (ncol(new.X) == 0) {
      polr.out <- polr(ordered(Y) ~ 1)
    }
    else {
      polr.out <- polr(ordered(Y) ~ new.X)
    }
    
    ## prior for transition matrix
    A0 <- trans.mat.prior(m=m, n=N, a=a, b=b)
    
    ## prior for beta error checking
    if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,K,1)  
    }
    if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
      cat("N(b0,B0) prior b0 not conformable.\n")
      stop("Please respecify and call MCMCoprobitChange() again.\n") 
    }   
    if(is.null(dim(B0))) {
      B0 <- B0 * diag(K)    
    }
    if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
      cat("N(b0,B0) prior B0 not conformable.\n")
      stop("Please respecify and call MCMCoprobitChange() again.\n")
    }
    marginal.likelihood  <- match.arg(marginal.likelihood)
    B0.eigenvalues <- eigen(B0)$values
     if (isTRUE(all.equal(min(B0.eigenvalues), 0))){
      if (marginal.likelihood != "none"){
        warning("Cannot calculate marginal likelihood with improper prior\n")
        marginal.likelihood <- "none"
      }
    }
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }
    
    ## to save time
    B0inv <- solve(B0)    
    gamma.start <- matrix(NA, ncat + 1, 1)
    gamma.start[1] <- -300
    gamma.start[2] <- 0
    gamma.start[3:ncat] <- (polr.out$zeta[2:(ncat - 1)] - polr.out$zeta[1]) * 0.588
    gamma.start[ncat + 1] <- 300
    
    ## initial values
    mle <- polr(Y ~ X[,-1])
    beta <- matrix(rep(c(mle$zeta[1], coef(mle)), ns), ns, , byrow=TRUE)
    ols <- lm(as.double(Y) ~ X-1)
    betalinearstart <- matrix(rep(coef(ols), ns), ns, , byrow=TRUE)
    P <-  trans.mat.prior(m=m, n=N, a=0.9, b=0.1)
    Sigmastart <- summary(ols)$sigma
    if (gamma.fixed==1){
      gamma <- gamma.start
      gamma.storage <-rep(0.0, nstore*gk)
    }
    else {
      gamma <- matrix(rep(gamma.start, ns), ns, ,  byrow=T)
      gamma.storage <- rep(0.0, nstore*ns*gk)
    }
    
    ## call C++ code to draw sample
    posterior <- .C("MCMCoprobitChange", 
                    betaout = as.double(rep(0.0, nstore*ns*K)), 
                    betalinearout = as.double(rep(0.0, nstore*ns*K)), 
                    gammaout = as.double(gamma.storage), 
                    Pout = as.double(rep(0.0, nstore*ns*ns)), 
                    psout = as.double(rep(0.0, N*ns)),
                    sout = as.double(rep(0.0, nstore*N)), 

                    Ydata = as.double(Y),
                    Xdata = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),

                    m = as.integer(m), 
                    ncat = as.integer(ncat),
                    
                    burnin = as.integer(burnin),           
                    mcmc = as.integer(mcmc), 
                    thin = as.integer(thin),
                    verbose = as.integer(verbose),

                    tunedata = as.double(tune),
                    lecuyer=as.integer(lecuyer), 
                    seedarray=as.integer(seed.array),
                    lecuyerstream=as.integer(lecuyer.stream),

                    betastart = as.double(beta),
                    betalinearstart = as.double(betalinearstart),
                    gammastart = as.double(gamma),
                    Pstart = as.double(P),    
                    sigmastart = as.double(Sigmastart),
                    
                    a = as.double(a),
                    b = as.double(b),
                    b0data = as.double(b0),
                    B0data = as.double(B0), 
                    A0data = as.double(A0), 
                    logmarglikeholder = as.double(0.0),
                    loglikeholder = as.double(0.0),
                    chib = as.integer(chib),
                    gammafixed= as.integer(gamma.fixed))                
               
    ## get marginal likelihood if Chib95
    if (chib==1){
      logmarglike <- posterior$logmarglikeholder
      loglike <- posterior$loglikeholder
    }
    else{
      logmarglike <- loglike <- 0
    }
    
    ## pull together matrix and build MCMC object to return
    beta.holder <- mcmc(matrix(posterior$betaout, nstore, ns*K))
    if (gamma.fixed==1){
      gamma.holder <- mcmc(matrix(posterior$gammaout, nstore, gk))
    }
    else {
      gamma.holder <- mcmc(matrix(posterior$gammaout, nstore, ns*gk))
    }
    P.holder    <- matrix(posterior$Pout, nstore, )
    s.holder    <- matrix(posterior$sout, nstore, )
    ps.holder   <- matrix(posterior$psout, N, )
    
    varnames(beta.holder)  <- sapply(c(1:ns),
                                     function(i){
                                       paste(c(xnames), "_regime", i, sep = "")
                                     })
    ## betalinear
    betalinear.holder <- mcmc(matrix(posterior$betalinearout, nstore, ns*K))
    varnames(betalinear.holder)  <- sapply(c(1:ns),
                                           function(i){
                                             paste(c(xnames), "_regime", i, sep = "")
                                           })
    gamma.holder <- gamma.holder[, as.vector(sapply(1:ns, function(i){gk*(i-1) + (3:(gk-1))}))]
    gamma.names <- paste("gamma", 3:(gk-1), sep="")
    varnames(gamma.holder)  <- sapply(c(1:ns),
                                      function(i){
                                        paste(gamma.names, "_regime", i, sep = "")
                                      })
    
    output <- mcmc(cbind(beta.holder, gamma.holder))
    attr(output, "title") <- "MCMCoprobitChange Posterior Sample"
    ## attr(output, "betalinear") <- mcmc(betalinear.holder)
    attr(output, "formula") <- formula
    attr(output, "y")       <- Y
    attr(output, "X")       <- X
    attr(output, "m")       <- m
    attr(output, "call")    <- cl
    attr(output, "logmarglike") <- logmarglike
    attr(output, "loglike") <- loglike
    attr(output, "prob.state") <- ps.holder/nstore
    attr(output, "s.store") <- s.holder
    return(output)
    
  }## end of MCMC function
