##########################################################################
## MCMCSVDreg.R samples from the posterior distribution of a Gaussian
## linear regression model in which the X matrix has been decomposed
## with an SVD. Useful for prediction when number of columns of X
## is (possibly much) greater than the number of rows of X.
##
## See West, Mike. 2000. "Bayesian Regression in the 'Large p, Small n'
##      Paradigm." Duke ISDS Discussion Paper 2000-22.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 9/9/2005
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

parse.formula.SVDreg <- function(formula, data, intercept){
  
  ## extract Y, X, and variable names for model formula and frame
  mt <- terms(formula, data=data)
  if(missing(data)) data <- sys.frame(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- NULL
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass ## for specialty handling of missing data
  
  mf[[1]] <- as.name("model.frame")

   mf <- eval(mf, sys.frame(sys.parent()))
  if (!intercept){
    attributes(mt)$intercept <- 0
  }

  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
  X <- as.matrix(X)         # X matrix
  Y <- as.matrix(model.response(mf, "numeric")) # Y matrix


  ## delete obs that are missing in X but potentially keep obs that are
  ## missing Y
  ## These are used to for the full SVD of X'
  keep.indic <- apply(is.na(X), 1, sum) == 0
  Y.full <- as.matrix(Y[keep.indic,])
  X.full <- X[keep.indic,]
  xvars.full <- dimnames(X.full)[[2]] # X variable names
  xobs.full  <- dimnames(X.full)[[1]] # X observation names
    
  return(list(Y.full, X.full, xvars.full, xobs.full))
  
}


"MCMCSVDreg" <-
  function(formula, data=NULL, burnin = 1000, mcmc = 10000,
           thin=1, verbose = 0, seed = NA, tau2.start = 1,
           g0 = 0, a0 = 0.001, b0 = 0.001, c0=2, d0=2, w0=1,
           beta.samp=FALSE, intercept=TRUE, ...) {
    
    # checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    
    # seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    ## form response and model matrices
    holder <- parse.formula.SVDreg(formula, data, intercept)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    obsnames <- holder[[4]]
    K <- ncol(X)  # number of covariates in unidentified model
    N <- nrow(X)  # number of obs (including obs for which predictions
                  # are required) N is also the length of gamma 
    Y.miss.indic <- as.numeric(is.na(Y))
    n.miss.Y <- sum(Y.miss.indic)
    Y[is.na(Y)] <- mean(Y, na.rm=TRUE)
    
    
    ## create SVD representation of t(X)
    svd.out <- svd(t(X)) # t(X) = A %*% D %*% F
    A <- svd.out$u
    D <- diag(svd.out$d)
    F <- t(svd.out$v)
    
    ## starting values and priors
    if (length(tau2.start) < N){
      tau2.start <- rep(tau2.start, length.out=N)
    }
    tau2.start <- matrix(tau2.start, N, 1)
    mvn.prior <- form.mvn.prior(g0, 0, N)
    g0 <- mvn.prior[[1]]
    c0 <- rep(c0, length.out=N)
    d0 <- rep(d0, length.out=N)
    check.ig.prior(a0, b0)
    for (i in 1:N){
      check.ig.prior(c0[i], d0[i])
    }
    w0 <- rep(w0, length.out=N)
    if (min(w0) < 0 | max(w0) > 1){
      cat("Element(s) of w0 not in [0, 1].\n")
      stop("Please respecify and call MCMCSVDreg again.\n")      
    }
    
    ## define holder for posterior sample
    if (beta.samp){
      sample <- matrix(data=0, mcmc/thin, n.miss.Y + 2*N + 1 + K)
    }
    else{
      sample <- matrix(data=0, mcmc/thin, n.miss.Y + 2*N + 1)
    }

    
    ## call C++ code to draw sample
    posterior <- .C("MCMCSVDreg",
                    sampledata = as.double(sample),
                    samplerow = as.integer(nrow(sample)),
                    samplecol = as.integer(ncol(sample)),
                    Y = as.double(Y),
                    Yrow = as.integer(nrow(Y)),
                    Ycol = as.integer(ncol(Y)),
                    Ymiss = as.integer(Y.miss.indic),
                    A = as.double(A),
                    Arow = as.integer(nrow(A)),
                    Acol = as.integer(ncol(A)),
                    D = as.double(D),
                    Drow = as.integer(nrow(D)),
                    Dcol = as.integer(ncol(D)),
                    F = as.double(F),
                    Frow = as.integer(nrow(F)),
                    Fcol = as.integer(ncol(F)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    tau2.start = as.double(tau2.start),
                    tau2row = as.integer(nrow(tau2.start)),
                    tau2col = as.integer(ncol(tau2.start)),
                    g0 = as.double(g0),
                    g0row = as.integer(nrow(g0)),
                    g0col = as.integer(ncol(g0)),
                    a0 = as.double(a0),
                    b0 = as.double(b0),
                    c0 = as.double(c0),
                    d0 = as.double(d0),
                    w0 = as.double(w0),
                    betasamp = as.integer(beta.samp),
                    PACKAGE="MCMCpack"
                    )
    
    ## pull together matrix and build MCMC object to return
    Y.miss.names <- NULL
    if (sum(Y.miss.indic) > 0){
      Y.miss.names <- paste("y", obsnames[Y.miss.indic==1], sep=".")
    }
    gamma.names <- paste("gamma", 1:N, sep=".")
    tau2.names <- paste("tau^2", 1:N, sep=".")
    beta.names <- paste("beta", xnames, sep=".")
    if (beta.samp){
      output <- form.mcmc.object(posterior,
                                 names=c(Y.miss.names, gamma.names, tau2.names,
                                   "sigma^2", beta.names),
                                 title="MCMCSVDreg Posterior Sample")
    }
    else{
      output <- form.mcmc.object(posterior,
                                 names=c(Y.miss.names, gamma.names, tau2.names,
                                   "sigma^2"),
                                 title="MCMCSVDreg Posterior Sample")      
    }
    return(output)
  }
