####################################################################
## test group-level breaks from panel residuals 
##
## written by Jong Hee Park 03/2009
## modified and integrated with other codes by JHP 07/2011
######################################################################

"testpanelGroupBreak" <-
  function(subject.id, time.id, resid, m=1, 
           mcmc=1000, burnin=1000, thin=1, verbose=0, 
           b0, B0, c0, d0, a = NULL, b = NULL, 
           seed = NA, marginal.likelihood = c("none", "Chib95"), ...){
    ## beta.start and sigma2.start are not arguments. OLS estimates will be used!
    ## subject.id is a numeric list indicating the group number. It should start from 1.
    ## time.id is a numeric list indicating the time, starting from 1. 
    ## subject.offset is the obs number from which a new subject unit starts
    ## time.offset is the obs number from which a new time unit starts when we stack data by time.id
    
    cl <- match.call()
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    ns <- m + 1
    nobs <- length(subject.id)
    newY <- matrix(resid, nobs, 1)
    newX <- matrix(1, nobs, 1)
    K <-  1
    
    ## Sort Data based on time.id
    oldTSCS <- cbind(time.id, subject.id, newY, newX)
    newTSCS <- oldTSCS[order(oldTSCS[,1]),]
    newYT <- as.matrix(newTSCS[,3])
    newXT <- as.matrix(newTSCS[,4])   
    b0 <- as.matrix(b0)
    B0 <- as.matrix(B0)
 
    nstore <- mcmc/thin
    nsubj <- length(unique(subject.id))
    ## subject.groupinfo matrix
    if (unique(subject.id)[1] != 1){
      stop("subject.id should start 1!")
    }
    subject.offset <- c(0, which(diff(sort(subject.id))==1)[-nsubj])
    nsubject.vec <- rep(NA, nsubj)
    for (i in 1:nsubj){
      nsubject.vec[i] <- sum(subject.id==unique(subject.id)[i])
    }
    subject.groupinfo <- cbind(unique(subject.id), subject.offset, nsubject.vec)
    
    ## time.groupinfo
    if(unique(time.id)[1] != 1){
      time.id <- time.id - unique(time.id)[1] + 1
      cat("time.id does not start from 1. So it is modified by subtracting the first unit of time.")
    }
    ntime <- max(nsubject.vec)## maximum time length
    ntime.vec <- rep(NA, ntime)
    for (i in 1:ntime){
      ntime.vec[i] <- sum(time.id==unique(time.id)[i])
    }
    time.offset <- c(0, which(diff(sort(time.id))==1)[-ntime])
    time.groupinfo <- cbind(unique(time.id), time.offset, ntime.vec)

    ## prior inputs
    if (m > 0){
      P0 <- trans.mat.prior(m=m, n=ntime, a=a, b=b)
    }
    else {
      P0 <- matrix(1, 1, 1)
    }
    betadraws <- matrix(data=0, nstore, ns*K)
    sigmadraws <- matrix(data=0, nstore, ns)
    psdraws <- matrix(data=0, ntime, ns)
    ols <- lm(newY ~ newX - 1)
    beta.start  <- rep(coef(ols)[1], ns)
    sigma2.start <- summary(ols)$sigma^2   

    ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    
    ## following MCMCregress, set chib as binary
    logmarglike <- loglik <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    ## call C++ code to draw sample
    posterior <- .C("HMMmultivariateGaussian",                    
                    betadata = as.double(betadraws),
                    betarow = as.integer(nrow(betadraws)),
                    betacol = as.integer(ncol(betadraws)),
                    
                    sigmadata = as.double(sigmadraws),
                    psout = as.double(psdraws),
                    nsubj = as.integer(nsubj), ntime = as.integer(ntime),
                    m = as.integer(m),
                    nobs = as.integer(nobs),
                    subjectid = as.integer(subject.id),
                    timeid = as.integer(time.id),
                     
                    Ydata = as.double(newY), Yrow = as.integer(nrow(newY)), Ycol = as.integer(ncol(newY)),
                    Xdata = as.double(newX), Xrow = as.integer(nrow(newX)), Xcol = as.integer(ncol(newX)),
                    YTdata = as.double(newYT), XTdata = as.double(newXT), 
                    burnin = as.integer(burnin), mcmc = as.integer(mcmc),
                    thin = as.integer(thin), verbose = as.integer(verbose),
                    
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                  
                    betastartdata = as.double(beta.start),
                    sigma2start = as.double(sigma2.start),

                    b0data = as.double(b0), B0data = as.double(B0),
                    c0 = as.double(c0), d0 = as.double(d0),                   
                    P0data = as.double(P0),
                    P0row = as.integer(nrow(P0)),
                    P0col = as.integer(ncol(P0)),                    

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
    xnames <-  sapply(c(1:K), function(i){paste("beta", i, sep = "")})
    output1 <- mcmc(data=beta.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
    output2 <- mcmc(data=sigma.samp, start=burnin+1, end=burnin + mcmc, thin=thin)
    if(m>1){
      varnames(output1)  <- sapply(c(1:ns),
                                   function(i){
                                     paste(xnames, "_regime", i, sep = "")
                                   })
      varnames(output2)  <- sapply(c(1:ns),
                                   function(i){
                                     paste("sigma2_regime", i, sep = "")
                                   })
    }
    output <- as.mcmc(cbind(output1, output2))
    
    attr(output, "title") <- "testpanelGroupBreak Posterior Sample"
    attr(output, "call")   <- cl
    attr(output, "y")       <- resid[1:ntime]
    attr(output, "m")       <- m
    attr(output, "nsubj")   <- nsubj
    attr(output, "ntime")   <- ntime
    if(m>0){
      ps.holder   <- matrix(posterior$psout, ntime, ns)
      attr(output, "prob.state") <- ps.holder/nstore
    }
    attr(output, "logmarglike") <- posterior$logmarglikeholder
    attr(output, "loglike") <- posterior$loglikeholder

    ## report the results in a simple manner
    
 
    return(output)
  }
