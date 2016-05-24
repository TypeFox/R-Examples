####################################################################
## HMM Gaussian Panel Random Effects Model
## y_it = x_it'*beta + w_it*bi + e_it, 
## e_it ~ N(0, sigma2)
##
## bi ~ N(0, D) : random effects coefficient
## D ~ IW(r0, R0) : covariance matrix for multiple random effects
## beta ~ N(b0, B0) : fixed effect coefficient
## sigma2 ~ IG(c0/2, d0/2) : random error
##
## written by Jong Hee Park 03/2009
## modified and integrated with other codes on 09/2011
######################################################################

"HMMpanelRE" <-
  function(subject.id, time.id, y, X, W, m=1, 
           mcmc=1000, burnin=1000, thin=1, verbose=0, 
           b0=0, B0=0.001, c0 = 0.001, d0 = 0.001, r0, R0, a = NULL, b = NULL, 
           seed = NA, beta.start = NA, sigma2.start = NA, D.start= NA, P.start = NA, 
           marginal.likelihood = c("none", "Chib95"), ...){
       
    cl <- match.call()
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    ns <- m + 1
    Y <- as.matrix(y)
    X <- as.matrix(X)
    W <- as.matrix(W)
    formula <- Y ~ X-1
    ols <- lm(formula)

    N <-  nrow(Y)
    K <-  ncol(X)
    Q <-  ncol(W)
    
    nobs <- nrow(X)
 
    ## Sort Data based on time.id
    oldTSCS <- cbind(time.id, subject.id, y, X, W)
    newTSCS <- oldTSCS[order(oldTSCS[,1]),]
    YT <- as.matrix(newTSCS[,3])
    XT <- as.matrix(newTSCS[,4:(4+K-1)])
    WT <- as.matrix(newTSCS[,(4+K):(4+K+Q-1)])
      
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    R0 <- as.matrix(R0)
 
    nstore <- mcmc/thin
    nsubj <- length(unique(subject.id))
    if (unique(subject.id)[1] != 1){
      stop("subject.id should start 1!")
    }
    ## subject.offset is the obs number from which a new subject unit starts
    subject.offset <- c(0, which(diff(sort(subject.id))==1)[-nsubj])
    ## col1: subj ID, col2: offset (C indexing), col3: #time periods in each subject
    nsubject.vec <- rep(NA, nsubj)
    for (i in 1:nsubj){
      nsubject.vec[i] <- sum(subject.id==unique(subject.id)[i])
    }
    subject.groupinfo <- cbind(unique(subject.id), subject.offset, nsubject.vec)

    
    ## time.groupinfo
    ## col1: time ID, col2: offset (C indexing), col3: # subjects in each time
    if(unique(time.id)[1] != 1){
      time.id <- time.id - unique(time.id)[1] + 1
      cat("time.id does not start from 1. So it is modified by subtracting the first unit of time.")
    }
    ntime <- max(nsubject.vec)## maximum time length
    ntime.vec <- rep(NA, ntime)
    for (i in 1:ntime){
      ntime.vec[i] <- sum(time.id==unique(time.id)[i])
    }
    ## time.offset is the obs number from which a new time unit starts when we stack data by time.id
    time.offset <- c(0, which(diff(sort(time.id))==1)[-ntime])
    time.groupinfo <- cbind(unique(time.id), time.offset, ntime.vec)

    
    ## prior inputs
    if (m > 0){
      P0 <- trans.mat.prior(m=m, n=ntime, a=a, b=b)
      ## initial values
      Pstart  <-  check.P(P.start, m, a=a, b=b)
    }
    else {
      Pstart <- P0 <- matrix(1, 1, 1)
    }
    if (is.na(beta.start[1])) {
      betastart  <- coef(ols)
    }
    else{
      betastart <- beta.start
    }
    if (is.na(sigma2.start[1])) {
      sigma2start <- summary(ols)$sigma^2
    }
    else{
      sigma2start <- sigma2.start
    }
    
    betadraws <- matrix(data=0, nstore, ns*K)
    sigmadraws <- matrix(data=0, nstore, ns)
    Ddraws <- matrix(data=0, nstore, ns*Q*Q)
    psdraws <- matrix(data=0, ntime, ns)
    sdraws <- matrix(data=0, nstore, ntime*ns)
    
    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    
    ## following MCMCregress, set chib as binary
    logmarglike <- loglik <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }
    
    ## call C++ code to draw sample
    posterior <- .C("HMMpanelRE",                    
                    betadata = as.double(betadraws),
                    betarow = as.integer(nrow(betadraws)),
                    betacol = as.integer(ncol(betadraws)),                    
                    sigmadata = as.double(sigmadraws),
                    Ddata = as.double(Ddraws),
                    psout = as.double(psdraws),
                    sout = as.double(sdraws),
                    
                    nsubj = as.integer(nsubj), ntime = as.integer(ntime), m = as.integer(m),
                    nobs = as.integer(nobs),
                    subjectid = as.integer(subject.id),
                    timeid = as.integer(time.id),
                     
                    Ydata = as.double(Y), Yrow = as.integer(nrow(Y)), Ycol = as.integer(ncol(Y)),
                    Xdata = as.double(X), Xrow = as.integer(nrow(X)), Xcol = as.integer(ncol(X)),
                    Wdata = as.double(W), Wrow = as.integer(nrow(W)), Wcol = as.integer(ncol(W)),
                    YTdata = as.double(YT), XTdata = as.double(XT), WTdata = as.double(WT),
                    burnin = as.integer(burnin), mcmc = as.integer(mcmc),
                    thin = as.integer(thin), verbose = as.integer(verbose),
                    
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                  
                    betastartdata = as.double(betastart),
                    sigma2start = as.double(sigma2start),
                    Pstart = as.double(Pstart),
 
                    b0data = as.double(b0), B0data = as.double(B0),
                    c0 = as.double(c0), d0 = as.double(d0),                   
                    r0 = as.integer(r0), R0data = as.double(R0),           

                    subject_groupinfodata = as.double(subject.groupinfo),
                    time_groupinfodata = as.double(time.groupinfo),

                    logmarglikeholder = as.double(0), 
                    loglikeholder = as.double(0), 
                    chib = as.integer(chib),
                    PACKAGE="MCMCpack"
                    )
    
    ## pull together matrix and build MCMC object to return
    beta.samp <- matrix(posterior$betadata,
                        posterior$betarow,
                        posterior$betacol)
    ## stored by the order of (11, 12, 13, 21, 22, 23)
    sigma.samp <- matrix(posterior$sigmadata,
                         posterior$betarow,
                         ns)    
    D.samp <- matrix(posterior$Ddata,
                     posterior$betarow,
                     Q*Q*ns)
    xnames <-  sapply(c(1:K), function(i){paste("beta", i, sep = "")})
    Dnames <-  sapply(c(1:(Q*Q)), function(i){paste("D", i, sep = "")})
    if (m == 0){
      output <- as.mcmc(cbind(beta.samp, sigma.samp, D.samp))
      names <- c(xnames, "sigma2", Dnames)
      varnames(output) <- as.list(names)
    }
    else{
      output1 <- mcmc(data=beta.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
      output2 <- mcmc(data=sigma.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
      output3 <- mcmc(data=D.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
      varnames(output1)  <- sapply(c(1:ns),
                                 function(i){
                                   paste(xnames, "_regime", i, sep = "")
                                 })
      varnames(output2)  <- sapply(c(1:ns),
                                   function(i){
                                     paste("sigma2_regime", i, sep = "")
                                   })
      varnames(output3)  <- sapply(c(1:ns),
                                   function(i){
                                     paste(Dnames, "_regime", i, sep = "")
                                   })
      
      output <- as.mcmc(cbind(output1, output2, output3))
      ps.holder   <- matrix(posterior$psout, ntime, ns)
      s.holder    <- matrix(posterior$sout, nstore, )
    }

    attr(output, "title") <- "HMMpanelRE Posterior Sample"
    attr(output, "call")   <- cl
    attr(output, "y")       <- y[1:ntime]
    attr(output, "X")       <- X[1:ntime, ]
    attr(output, "m")       <- m
    attr(output, "nsubj")   <- nsubj
    attr(output, "ntime")   <- ntime
    if (m > 0){
      attr(output, "s.store") <- s.holder
      attr(output, "prob.state") <- ps.holder/nstore
    }
    attr(output, "logmarglike") <- posterior$logmarglikeholder
    attr(output, "loglike") <- posterior$loglikeholder 
    
    return(output)
  }
