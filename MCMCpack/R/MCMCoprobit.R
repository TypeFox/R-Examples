##########################################################################
## sample from the posterior distribution of an ordered probit model
## via the data augmentation approach of Cowles (1996)
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 1/25/2003
## Modified to meet new developer specification 7/26/2004 KQ
## Modified for new Scythe and rngs 7/26/2004 KQ
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################


"MCMCoprobit" <-
  function(formula, data = parent.frame(), burnin = 1000, mcmc = 10000,
           thin = 1, tune = NA, tdf = 1, verbose = 0, seed = NA, beta.start = NA,
           b0 = 0, B0 = 0, a0 = 0, A0 = 0, mcmc.method = c("Cowles", "AC"), ...) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    ## extract X, Y, and variable names from the model formula and frame 
    call <- match.call() 	  
    mt <- terms(formula, data=data)
    if(missing(data)) data <- sys.frame(sys.parent())
    mf <- match.call(expand.dots = FALSE)
    mf$burnin <- mf$mcmc <- mf$b0 <- mf$B0 <- mf$a0 <- mf$A0 <- NULL
    mf$thin <- mf$... <- mf$tune <- mf$tdf <- mf$verbose <- mf$seed <- NULL
    mf$beta.start <- mf$mcmc.method <- NULL
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    vars <- as.character(attr(mt, "variables"))[-1] # y varname and x varnames
    
    ## null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)# else NULL
    X.names <- dimnames(X)[[2]]
    Y    <- model.response(mf, "numeric")
    Y    <- factor(Y, ordered=TRUE)
    ncat <- nlevels(Y)             # number of categories of y
    cat  <- levels(Y)              # values of categories of y
    N <- nrow(X)	               # number of observations
    K <- ncol(X)	               # number of covariates
    if (length(Y) != N){
      cat("X and Y do not have same number of rows.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")      
    }
    
    ## convert data to matrices to be passed
    Y <- as.matrix(as.integer(Y))
    X <- as.matrix(X)
    
    ## check tuning parameter
    if (is.na(tune)){
      tune <- 0.05/ncat
    }
    
    xint <- match("(Intercept)", colnames(X), nomatch=0)
    if (xint > 0){
      new.X <- X[, -xint, drop=FALSE]
    }
    else warning("An intercept is needed and assumed in MCMCoprobit()\n.")
    if (ncol(new.X) == 0){
      polr.out <- polr(ordered(Y)~1)
    }
    else {
      polr.out <- polr(ordered(Y)~new.X)
    }
    
    ## starting values for beta error checking
    if (is.na(beta.start)){
      beta.start <- matrix(0, K, 1)
      beta.start[1] <- -.588 * polr.out$zeta[1]
      if( ncol(new.X) > 0){
        beta.start[2:K] <- .588 * coef(polr.out)
      }
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- beta.start * matrix(1,K,1)  
    }
    else if((dim(beta.start)[1] != K) || (dim(beta.start)[2] != 1)) {
      cat("Starting value for beta not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }
    
    ## prior for beta error checking
    if(is.null(dim(b0))) {
      b0 <- b0 * matrix(1,K,1)  
    }
    if((dim(b0)[1] != K) || (dim(b0)[2] != 1)) {
      cat("N(b0,B0) prior b0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n") 
    }   
    if(is.null(dim(B0))) {
      B0 <- B0 * diag(K)    
    }
    if((dim(B0)[1] != K) || (dim(B0)[2] != K)) {
      cat("N(b0,B0) prior B0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }

    ## prior for alpha error checking
    if(is.null(dim(a0))) {
      a0 <- a0 * matrix(1, ncat-1, 1)  
    }
    if((dim(a0)[1] != ncat-1) || (dim(a0)[2] != 1)) {
      cat("N(a0,A0) prior a0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n") 
    }   
    if(is.null(dim(A0))) {
      A0 <- A0 + diag(ncat - 1)    
    }
    if((dim(A0)[1] != ncat - 1) || (dim(A0)[2] != ncat - 1)) {
      cat("N(a0, A0) prior A0 not conformable.\n")
      stop("Please respecify and call MCMCoprobit() again.\n")
    }
    
    ## form gamma starting values (note: not changeable)
    gamma <- matrix(NA,ncat+1,1)
    gamma[1] <- -300
    gamma[2] <- 0
    gamma[3:ncat] <- (polr.out$zeta[2:(ncat-1)] - polr.out$zeta[1])*.588
    gamma[ncat+1] <- 300
    
    ## posterior sample
    sample <- matrix(data=0, mcmc/thin, K + ncat + 1)
    
    ## call C++ code to draw sample
    nY <- as.matrix(as.numeric(Y))
    
    ## mcmc.method
    cowles <- as.integer(1)
    if(mcmc.method[1]=="AC") {cowles <- as.integer(0)}
    
    ## form the tuning parameter
    tune <- vector.tune(tune, ncat-1)
    posterior <- NULL
    
    auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCoprobit",
                     sample.nonconst=sample, Y=as.integer(Y), nY=nY, X=X, 
                     burnin=as.integer(burnin),
                     mcmc=as.integer(mcmc), thin=as.integer(thin),
                     tune=tune, tdf=as.double(tdf), 
                     lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose), beta=beta.start,
                     gamma=gamma, b0=b0, B0=B0, a0=a0, A0=A0, 
                     cowles=as.integer(cowles)) 
    
    ## put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol, byrow=FALSE)
    if(mcmc.method[1]=="AC"){ 
      sample[ , 1] <- sample[, 1] - sample[, K+2] ## post-MCMC normalization 
      sample[ , (K+2):(K+ncat)] <- sample[ , (K+2):(K+ncat)] - sample[, K+2] ## post-MCMC normalization 
    }
    sample <- sample[, c(1:K, (K+3):(K+ncat))]
    
    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    xnames <- c(X.names, paste("gamma", 2:(ncat-1), sep=""))
    varnames(output) <- xnames
    attr(output, "title") <- "MCMCoprobit Posterior Sample"
    
    return(output)
  }
