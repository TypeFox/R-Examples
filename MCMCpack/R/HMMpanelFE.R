####################################################################
## HMM Gaussian Panel with Time Varying Intercepts  
##
## first written by Jong Hee Park on 02/2009
## revised for MCMCpack inclusion on 09/2011
######################################################################
"HMMpanelFE" <-
  function(subject.id, y, X, m,
           mcmc=1000, burnin=1000, thin=1, verbose=0, 
           b0=0, B0=0.001, c0 = 0.001, d0 = 0.001, delta0=0, Delta0=0.001,
           a = NULL, b = NULL, seed = NA, ...){
    ## m should be a vector with a number of breaks for each group
    ## id is a numeric list
    ## p is a lag order 
    ## offset is the first time period from which each group starts

    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## Data
    Y <- as.matrix(y);
    X <- as.matrix(cbind(1, X)) ## the intercept is not reported
    N <-  nrow(Y);
    K <-  ncol(X)
    
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]
    
    mvn.prior <- form.mvn.prior(delta0, Delta0, 1)
    delta0 <- mvn.prior[[1]]
    Delta0 <- mvn.prior[[2]]
    
    nstore <- mcmc/thin
    nsubj <- length(unique(subject.id))
    
    ## groupinfo matrix
    ## col1: subj ID, col2: offset (first time C indexing), col3: #time periods
    if (unique(subject.id)[1] != 1){
      stop("subject.id should start 1!")
    }
    subject.offset <- c(0, which(diff(sort(subject.id))==1)[-nsubj])
    ## col1: subj ID, col2: offset (C indexing), col3: #time periods in each subject
    nsubject.vec <- rep(NA, nsubj)
    for (i in 1:nsubj){
      nsubject.vec[i] <- sum(subject.id==unique(subject.id)[i])
    }
    subject.groupinfo <- cbind(unique(subject.id), subject.offset, nsubject.vec)
    
    ## maximum time length
    ntime <- max(nsubject.vec)
    m.max <- max(m)
    m.min <- min(m)
    
    ## prior inputs
    P0data <- NULL
    Pstart <- NULL
    for (i in 1:nsubj){
      if(m[i] == 0){
        P0current <- 1
        Pstartcurrent <- 1
      }
      else{
        P0current <- trans.mat.prior(m=m[i], n=nsubject.vec[i], a=a, b=b)
        Pstartcurrent <- trans.mat.prior(m=m[i], n=nsubject.vec[i], a=.9, b=.1)
      }
      P0data <- c(P0data, P0current)
      Pstart <- c(Pstart, Pstartcurrent) 
    }

    ## starting values
    ols <- lm(Y~X-1)
    beta.start <- coef(ols)
    sigma2.start <- summary(ols)$sigma^2
    deltastart  <- NULL
    Ytilde <- Y - X%*%beta.start
    deltaformula <- Ytilde ~ 1 ## without intercept
    for (i in 1:nsubj){
      deltacurrent <- rep(as.vector(coef(lm(Ytilde ~ 1))), m[i] + 1)
      deltastart <- c(deltastart, deltacurrent) 
    }

    ## Storage
    totalstates0 <- sum(m+1)
    betadraws0 <- matrix(0, nstore, K)
    deltadraws0 <- matrix(data=0, nstore, totalstates0)
    sigmadraws0 <- matrix(data=0, nstore, totalstates0)
    statedraws0 <- matrix(data=0, nstore, totalstates0)
    
    ## call C++ code to draw sample
    posterior <- .C("HMMpanelFE",
                    deltadraws = as.double(deltadraws0),
                    sigmadraws = as.double(sigmadraws0),
                    statedraws = as.double(statedraws0),
                     
                    betadraws = as.double(betadraws0),
                    betarow = as.integer(nrow(betadraws0)),
                    betacol = as.integer(ncol(betadraws0)),
                    totalstates = as.integer(totalstates0),             
                                 
                    nsubj = as.integer(nsubj),
                    ntime = as.integer(ntime),
                    nobs = as.integer(N),
                    subjectid = as.integer(subject.id),
                    
                    m = as.integer(m),
                    mmax = as.integer(m.max),
                    mmin = as.integer(m.min),

                    Ydata = as.double(Y),
                    Yrow = as.integer(nrow(Y)),
                    Ycol = as.integer(ncol(Y)),

                    Xdata = as.double(X),
                    Xrow = as.integer(nrow(X)),
                    Xcol = as.integer(ncol(X)),
                    
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    verbose = as.integer(verbose),
                    
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                  
                    betastartdata = as.double(beta.start),
                    sigma2start = as.double(sigma2.start),
                    deltastartdata = as.double(deltastart),
                    deltastartrow = as.integer(length(deltastart)),
                    
                    b0data = as.double(b0),                 
                    B0data = as.double(B0),

                    delta0data = as.double(delta0),
                    Delta0data = as.double(Delta0),
                    
                    c0 = as.double(c0),
                    d0 = as.double(d0),                   
 
                    P0data = as.double(P0data),
                    P0row = as.integer(length(P0data)),
                    Pstartdata = as.double(Pstart),

                    subject_groupinfodata = as.double(subject.groupinfo),
                    PACKAGE="MCMCpack")
    
    ## pull together matrix and build MCMC object to return
    betadraws <- matrix(posterior$betadraws,
                        posterior$betarow,
                        posterior$betacol)
    
    sigma.samp <- as.mcmc(matrix(posterior$sigmadraws, nstore, totalstates0))
    delta.samp <- as.mcmc(matrix(posterior$deltadraws, nstore, totalstates0))
    state.samp <- as.mcmc(matrix(posterior$statedraws, nstore, totalstates0))
    
    ## output <- mcmc(betadraws, start=burnin+1, end=burnin+mcmc, thin=thin)
    output <- as.mcmc(betadraws[, -1])## drop the intercept  
    attr(output, "title") <- "HMMpanelFE Sample"
    attr(output, "m")  <- m
    attr(output, "sigma") <- sigma.samp
    attr(output, "state") <- state.samp
    attr(output, "delta") <- delta.samp
    return(output)
  }
