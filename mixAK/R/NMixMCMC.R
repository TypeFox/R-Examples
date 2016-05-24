##
##  PURPOSE:   (Reversible jump) MCMC for a normal mixture model
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:                  29/10/2007
##
##  IMPORTANT MODIFICATIONS:  03/11/2008  (computation of penalized expected deviance added)
##                                        (basic support for parallel computation on multicore
##                                         CPUs added using packages snowfall/snow)
##                            08/02/2013  snow/snowfall support for parallel computation replaced by parallel package
##                            26/03/2015  dependence of weights on a categorical covariate started
##                                        to be implemented
##
##  MINOR MODIFICATION:       26/04/2009  computation of initial values of censored observations
##                                        and init2$mu changed (variable tmpsd computed in other way
##                                        as suggested by referee of the paper)
##
##  FUNCTIONS:  NMixMCMC
##
## ======================================================================

## *************************************************************
## NMixMCMC
## *************************************************************
NMixMCMC <- function(y0, y1, censor, x_w, scale, prior,
                     init, init2, RJMCMC,
                     nMCMC = c(burn = 10, keep = 10, thin = 1, info = 10),
                     PED, keep.chains = TRUE, onlyInit = FALSE, dens.zero = 1e-300, parallel = FALSE, cltype)
{
  thispackage <- "mixAK"

  EMin <- -5           ## exp(-(D1+D2)) = exp(-EMin) when computing importance sampling weights
                       ## if D1 + D2 < EMin, where D1 = log(f(y|theta1)), D2 = log(f(y|theta2))
    ## -> these constants are passed to .C("NMix_PED")
    ## -> dens.zero is also used in NMixMCMCwrapper to compute D.bar
   
  
  ########## ========== Data ========== ##########
  ########## ========================== ##########
  dd <- NMixMCMCdata(y0 = y0, y1 = y1, censor = censor)
  rm(list=c("y0", "y1", "censor"))
    ### use dd$y0, dd$y1, dd$censor instead  

  ######### ========= Covariates for mixture weights (if any) =========== ############
  if (missing(x_w)){
    dd$x_w <- rep(0, dd$n)
    dd$fx_w <- factor(dd$x_w)
    dd$nx_w <- 1
    dd$lx_w <- levels(dd$fx_w)
  }else{
    if (length(x_w) != dd$n) stop("argument x_w should be a vector of length ", dd$n, " (it has a length ", length(x_w),").")
    if (any(is.na(x_w))) stop("NA's in a vector x_w are not allowed.")
    if (is.factor(x_w)) dd$fx_w <- x_w else dd$fx_w <- factor(x_w)
    dd$x_w <- as.numeric(dd$fx_w) - 1   ## values 0, 1, ...
    dd$lx_w <- levels(dd$fx_w)
    dd$nx_w <- length(dd$lx_w)
    rm(list = "x_w")
  }    

  
  ######### =========== Temporar initial values (equal to lower limit of right-censored, =========== #########
  ######### =========== upper limit of left-censored                                     =========== #########
  ######### =========== and midpoint of interval-censored observation)                   =========== #########
  ########## ======================================================================================= #########
  tmpinity <- dd$y0
  if (dd$are.Interval) tmpinity[dd$censor == 3] <- (dd$y0[dd$censor == 3] + dd$y1[dd$censor == 3])/2      
  
  
  ########## ========== Prior and Initial Values (chain 1), Scaling the data ========== ##########
  ########## ========================================================================== ##########
  if (missing(prior)) stop("prior must be given")
  if (!is.list(prior)) stop("prior must be a list")
  if (!length(prior)) stop("prior has a zero length")

  if (missing(init)) init <- list()
  if (!is.list(init)) stop("init must be a list")
  
  inprior <- names(prior)
  ipriorK   <- match("priorK",   inprior, nomatch=NA)
  ipriormuQ <- match("priormuQ", inprior, nomatch=NA)
  iKmax     <- match("Kmax",     inprior, nomatch=NA)
  ilambda   <- match("lambda",   inprior, nomatch=NA)
  idelta    <- match("delta",    inprior, nomatch=NA)
  ixi       <- match("xi",       inprior, nomatch=NA)  
  ice       <- match("ce",       inprior, nomatch=NA)
  iD        <- match("D",        inprior, nomatch=NA)
  izeta     <- match("zeta",     inprior, nomatch=NA)
  ig        <- match("g",        inprior, nomatch=NA)
  ih        <- match("h",        inprior, nomatch=NA)

  ininit    <- names(init)
  iy        <- match("y",        ininit, nomatch=NA)    
  iK        <- match("K",        ininit, nomatch=NA)
  iw        <- match("w",        ininit, nomatch=NA)
  imu       <- match("mu",       ininit, nomatch=NA)
  iSigma    <- match("Sigma",    ininit, nomatch=NA)
  iLi       <- match("Li",       ininit, nomatch=NA)  
  igammaInv <- match("gammaInv", ininit, nomatch=NA)
  ir        <- match("r",        ininit, nomatch=NA)    
  
  ##### integer prior:  Kprior
  ##### ----------------------------------------------------
  if (is.na(ipriorK)) prior$priorK <- "fixed"
  if (length(prior$priorK) != 1) stop("prior$priorK must be of length 1")
  CpriorK <- pmatch(prior$priorK, table=c("fixed", "uniform", "tpoisson"), nomatch=0) - 1
  if (CpriorK == -1) stop("prior$priorK must be one of fixed/uniform/tpoisson")

  if (prior$priorK != "fixed" & dd$nx_w > 1) stop("covariates on mixture weights not allowed if K is not fixed.")
  
  ##### PED
  ##### ----------------------------------------------------
  if (missing(PED)){
    if (prior$priorK == "fixed") PED <- TRUE
    else                         PED <- FALSE
  }  

  
  ##### integer prior:  priormuQ
  ##### ----------------------------------------------------  
  if (is.na(ipriormuQ)) prior$priormuQ <- "independentC"
  if (length(prior$priormuQ) != 1) stop("prior$priormuQ must be of length 1")  
  CpriormuQ <- pmatch(prior$priormuQ, table=c("naturalC", "independentC"), nomatch=0) - 1
  if (CpriormuQ == -1) stop("prior$priormuQ must be one of naturalC/independentC")

  ##### integer prior:  Kmax
  ##### ----------------------------------------------------  
  if (is.na(iKmax)) stop("prior$Kmax must be given")
  if (length(prior$Kmax) != 1) stop("prior$Kmax must be of length 1")
  if (is.na(prior$Kmax)) stop("NA in prior$Kmax")
  if (prior$Kmax <= 0) stop("prior$Kmax must be positive")
  CKmax <- as.numeric(prior$Kmax)
    
  ##### double prior:  lambda
  ##### ----------------------------------------------------  
  if (CpriorK == 2){    ## truncated Poisson prior for K
    if (is.na(ilambda)) stop("prior$lambda must be given when prior$priorK = tpoisson")
    if (length(prior$lambda) != 1) stop("prior$lambda must be of length 1")
    if (is.na(prior$lambda)) stop("NA in prior$lambda")
    if (prior$lambda <= 0) stop("prior$lambda must be positive")    
  }else{
    prior$lambda <- 0
  }
  Clambda <- as.numeric(prior$lambda)
  names(Clambda) <- "lambda"
  
  ##### double prior:  delta
  ##### ----------------------------------------------------  
  if (is.na(idelta)) prior$delta <- 1
  if (length(prior$delta) != 1) stop("prior$delta must be of length 1")
  if (is.na(prior$delta)) stop("NA in prior$delta")  
  if (prior$delta <= 0) stop("prior$delta must be positive")
  Cdelta <- as.numeric(prior$delta)
  names(Cdelta) <- "delta"

  
  ##### init:  y
  ##### ----------------------------------------------------
  if (is.na(iy)){
    init$y <- NMixMCMCinity(y0=dd$y0, y1=dd$y1, censor=dd$censor, sd.init=sd(tmpinity),
                            are.Censored=dd$are.Censored, are.Right=dd$are.Right, are.Exact=dd$are.Exact, are.Left=dd$are.Left, are.Interval=dd$are.Interval,
                            p=dd$p, n=dd$n, random=FALSE)    
  }else{
    init$y <- NMixMCMCinity(y0=dd$y0, y1=dd$y1, censor=dd$censor, sd.init=sd(tmpinity),
                            are.Censored=dd$are.Censored, are.Right=dd$are.Right, are.Exact=dd$are.Exact, are.Left=dd$are.Left, are.Interval=dd$are.Interval,
                            p=dd$p, n=dd$n, inity=init$y)    
  }  
  
    
  ##### scale:  Scaling the data
  ##### ------------------------------------------------------
  if (missing(scale)){                   ### REMARK: dd$y0, dd$y1 and init$y are already of class matrix (even if p = 1)
    SHIFT <- apply(init$y, 2, mean)
    SCALE <- apply(init$y, 2, sd)
    scale <- list(shift=SHIFT, scale=SCALE)
    rm(list=c("SHIFT", "SCALE"))
  }  
  if (!is.list(scale)) stop("scale must be a list")
  if (length(scale) != 2) stop("scale must have 2 components")
  inscale <- names(scale)  
  iscale.shift <- match("shift", inscale, nomatch=NA)
  iscale.scale <- match("scale", inscale, nomatch=NA)
  if (is.na(iscale.shift)) stop("scale$shift is missing")
  if (length(scale$shift) == 1) scale$shift <- rep(scale$shift, dd$p)
  if (length(scale$shift) != dd$p) stop(paste("scale$shift must be a vector of length ", dd$p, sep=""))    
  if (is.na(iscale.scale)) stop("scale$scale is missing")
  if (length(scale$scale) == 1) scale$scale <- rep(scale$scale, dd$p)
  if (length(scale$scale) != dd$p) stop(paste("scale$scale must be a vector of length ", dd$p, sep=""))
  if (any(scale$scale <= 0)) stop("all elements of scale$scale must be positive")

  z0    <- (dd$y0 - matrix(rep(scale$shift, dd$n), ncol=dd$p, byrow=TRUE))/matrix(rep(scale$scale, dd$n), ncol=dd$p, byrow=TRUE)
  z1    <- (dd$y1 - matrix(rep(scale$shift, dd$n), ncol=dd$p, byrow=TRUE))/matrix(rep(scale$scale, dd$n), ncol=dd$p, byrow=TRUE)    
  tmpinitz <- (tmpinity - matrix(rep(scale$shift, dd$n), ncol=dd$p, byrow=TRUE))/matrix(rep(scale$scale, dd$n), ncol=dd$p, byrow=TRUE)
  
  if (dd$p == 1){
    initz <- (init$y - scale$shift)/scale$scale
    zBar <- mean(initz)
    zMin <- min(initz)
    zMax <- max(initz)
  }else{
    initz <- (init$y - matrix(rep(scale$shift, dd$n), ncol=dd$p, byrow=TRUE))/matrix(rep(scale$scale, dd$n), ncol=dd$p, byrow=TRUE)
    zBar <- apply(initz, 2, mean)    ### will be used to determine prior$xi if this not given by the user
    zMin <- apply(initz, 2, min)     
    zMax <- apply(initz, 2, max)
  }
  zBar[abs(zBar) < 1e-14] <- 0  
  zVar <- var(initz)                 ### will be used to determine init$Sigma and init$gammaInv if these not given by the user
  zVar[abs(zVar - 1) < 1e-14] <- 1  
  zR <- zMax - zMin                  ### will be used to determine prior$D if this not given by the user
  zMid <- 0.5*(zMin + zMax)    
  
  
  ##### double prior:  xi
  ##### ----------------------------------------------------  
  if (is.na(ixi)) prior$xi <- matrix(rep(zMid, CKmax), nrow=CKmax, ncol=dd$p, byrow=TRUE)
  if (any(is.na(prior$xi))) stop("NA in prior$xi")  
  if (dd$p == 1){
    if (length(prior$xi) == 1) prior$xi <- rep(prior$xi, CKmax)                                                ## common prior$xi for all mixture components
    if (length(prior$xi) != CKmax) stop(paste("prior$xi must be of length ", CKmax, sep=""))
    prior$xi <- as.numeric(prior$xi)
    names(prior$xi) <- paste("xi", 1:CKmax, sep="")    
    Cxi <- prior$xi
  }else{
    if (length(prior$xi) == dd$p) prior$xi <- matrix(rep(as.numeric(prior$xi), each=CKmax), nrow=CKmax, ncol=dd$p)   ## common prior$xi for all mixture components
    if (CKmax == 1) prior$xi <- matrix(as.numeric(prior$xi), nrow=1)    
    if (!is.matrix(prior$xi)) stop("prior$xi must be a matrix")
    if (ncol(prior$xi) != dd$p) stop(paste("prior$xi must have ", dd$p, " columns", sep=""))
    if (nrow(prior$xi) != CKmax) stop(paste("prior$xi must have ", CKmax, " rows", sep=""))
    rownames(prior$xi) <- paste("j", 1:CKmax, sep="")
    colnames(prior$xi) <- paste("m", 1:dd$p, sep="")
    Cxi <- as.numeric(t(prior$xi))
    names(Cxi) <- paste("xi", rep(1:CKmax, each=dd$p), ".", rep(1:dd$p, CKmax), sep="")
  }  
  if (any(is.na(Cxi))) stop("NA in prior$xi")

  ##### double prior:  ce
  ##### ----------------------------------------------------  
  if (CpriormuQ == 0){    ## natural conjugate prior for (mu, Q)
    if (is.na(ice)) prior$ce <- rep(1, CKmax)
    if (length(prior$ce) == 1) prior$ce <- rep(prior$ce, CKmax)
    if (length(prior$ce) != CKmax) stop(paste("prior$ce must be of length ", CKmax, sep=""))
    if (any(is.na(prior$ce))) stop("NA in prior$ce")
    if (any(prior$ce <= 0)) stop("prior$ce must be positive")
    prior$ce <- as.numeric(prior$ce)
  }else{
    prior$ce <- rep(0, CKmax)
  }  
  Cce <- prior$ce
  names(Cce) <- names(prior$ce) <- paste("c", 1:CKmax, sep="")

  ##### double prior:  D
  ##### ----------------------------------------------------  
  if (CpriormuQ == 1){    ## independent conjugate prior for (mu, Q)
    if (is.na(iD)){
      if (dd$p == 1) prior$D <- rep(zR^2, CKmax)
      else           prior$D <- t(matrix(rep(diag(zR^2), CKmax), nrow=dd$p, ncol=CKmax*dd$p))
    }
    if (any(is.na(prior$D))) stop("NA in prior$D")    
    if (dd$p == 1){
      if (length(prior$D) == 1) prior$D <- rep(prior$D, CKmax)
      if (length(prior$D) != CKmax) stop(paste("prior$D must be of length ", CKmax, sep=""))
      prior$D <- as.numeric(prior$D)
      names(prior$D) <- paste("D", 1:CKmax, sep="")      
      if (any(prior$D <= 0)) stop("prior$D must be positive")
      CDinv <- 1/prior$D
      names(CDinv) <- paste("Dinv", 1:CKmax, sep="")      
    }else{
      if (!is.matrix(prior$D)) stop("prior$D must be a matrix")
      if (ncol(prior$D) != dd$p) stop(paste("prior$D must have ", dd$p, " columns", sep=""))
      if (nrow(prior$D) == dd$p){
        if (any(prior$D[lower.tri(prior$D)] != t(prior$D)[lower.tri(prior$D)])) stop("prior$D must be a symmetric matrix")
        err <- try(Dinv <- chol(prior$D), silent=TRUE)
        if (class(err) == "try-error") stop("Cholesky decomposition of prior$D failed")
        Dinv <- chol2inv(Dinv)
        CDinv <- rep(Dinv[lower.tri(Dinv, diag=TRUE)], CKmax)        
        prior$D <- matrix(rep(as.numeric(t(prior$D)), CKmax), nrow=dd$p*CKmax, ncol=dd$p, byrow=TRUE)
      }else{  
        if (nrow(prior$D) != CKmax*dd$p) stop(paste("prior$D must have ", CKmax, " times ", dd$p, " rows", sep=""))
        CDinv <- numeric(0)
        for (j in 1:CKmax){
          Dinv <- prior$D[((j-1)*dd$p+1):(j*dd$p),]
          if (any(Dinv[lower.tri(Dinv)] != t(Dinv)[lower.tri(Dinv)])) stop(paste(j, "-th block of prior$D is not symmetric", sep=""))
          err <- try(Dinv <- chol(Dinv), silent=TRUE)
          if (class(err) == "try-error") stop(paste("Cholesky decomposition of the ", j, "-th block of prior$D failed", sep=""))
          Dinv <- chol2inv(Dinv)
          CDinv <- c(CDinv, Dinv[lower.tri(Dinv, diag=TRUE)])
        }  
      }
      colnames(prior$D) <- paste("m", 1:dd$p, sep="")
      rownames(prior$D) <- paste("j", rep(1:CKmax, each=dd$p), ".", rep(1:dd$p, CKmax), sep="")
      names(CDinv) <- paste("Dinv", rep(1:CKmax, each=dd$LTp), rep(dd$naamLTp, CKmax), sep="")
    }  
  }else{
    if (dd$p == 1){
      prior$D <- rep(1, CKmax)
      names(prior$D) <- paste("D", 1:CKmax, sep="")
      CDinv <- 1/prior$D
      names(CDinv) <- paste("Dinv", 1:CKmax, sep="")      
    }else{
      prior$D <- matrix(rep(as.numeric(diag(dd$p)), CKmax), nrow=dd$p*CKmax, ncol=dd$p, byrow=TRUE)
      colnames(prior$D) <- paste("m", 1:dd$p, sep="")
      rownames(prior$D) <- paste("j", rep(1:CKmax, each=dd$p), ".", rep(1:dd$p, CKmax), sep="")

      Dinv <- diag(dd$p)
      CDinv <- rep(Dinv[lower.tri(Dinv, diag=TRUE)], CKmax)
      names(CDinv) <- paste("Dinv", rep(1:CKmax, each=dd$LTp), rep(dd$naamLTp, CKmax), sep="")      
    }  
  }  

  ##### double prior:  zeta
  ##### ----------------------------------------------------  
  if (is.na(izeta)) prior$zeta <- dd$p + 1
  if (length(prior$zeta) != 1) stop("prior$zeta must be of length 1")  
  if (is.na(prior$zeta)) stop("NA in prior$zeta")
  if (prior$zeta <= dd$p - 1) stop(paste("prior$zeta must be higher than ", dd$p - 1, sep=""))
  Czeta <- as.numeric(prior$zeta)
  names(Czeta) <- "zeta"
  
  ##### double prior:  g
  ##### ----------------------------------------------------  
  if (is.na(ig)) prior$g <- rep(0.2, dd$p)
  if (length(prior$g) == 1) prior$g <- rep(prior$g, dd$p)
  if (length(prior$g) != dd$p) stop(paste("prior$g must be of length ", dd$p, sep=""))  
  if (any(is.na(prior$g))) stop("NA in prior$g")
  if (any(prior$g <= 0)) stop("prior$g must be positive")
  Cg <- as.numeric(prior$g)
  names(Cg) <- paste("g", 1:dd$p, sep="")
  
  ##### double prior:  h
  ##### ----------------------------------------------------
  if (is.na(ih)) prior$h <- 10/(zR^2)
  if (length(prior$h) == 1) prior$h <- rep(prior$h, dd$p)
  if (length(prior$h) != dd$p) stop(paste("prior$h must be of length ", dd$p, sep=""))  
  if (any(is.na(prior$h))) stop("NA in prior$h")
  if (any(prior$h <= 0)) stop("prior$h must be positive")
  Ch <- as.numeric(prior$h)
  names(Ch) <- paste("h", 1:dd$p, sep="")

  
  ##### integer, double prior:  concetenate
  ##### ----------------------------------------------------
  Cinteger <- c(CpriorK, CpriormuQ, CKmax)
  names(Cinteger) <- c("priorK", "priormuQ", "Kmax")  
  Cdouble <- c(Clambda, Cdelta, Cxi, Cce, CDinv, Czeta, Cg, Ch)

  
  ##### init:  K
  ##### ----------------------------------------------------  
  if (is.na(iK)){
    if (prior$priorK == "fixed") init$K <- CKmax
    else                         init$K <- 1
  }
  if (prior$priorK == "fixed") init$K <- CKmax  
  if (length(init$K) != 1) stop("init$K must be of length 1")
  if (is.na(init$K)) stop("NA in init$K")
  if (init$K <= 0 | init$K > CKmax) stop("init$K out of the range")

  
  ##### init:  w
  ##### ----------------------------------------------------
  if (is.na(iw)){
    init$w <- rep(rep(1, init$K)/init$K, dd$nx_w)    
  }  
  init$w <- as.numeric(init$w)
  if (length(init$w) == CKmax * dd$nx_w & CKmax > init$K) init$w <- init$w[1:(init$K * dd$nx_w)]
  if (dd$nx_w > 1) names(init$w) <- paste("w", rep(1:init$K, dd$nx_w), ".", rep(1:dd$nx_w, each = init$K), sep="") else names(init$w) <- paste("w", 1:init$K, sep="")
  if (any(is.na(init$w))) stop("NA in init$w")  
  if (length(init$w) != init$K * dd$nx_w) stop(paste("init$w must be of length ", init$K * dd$nx_w, sep=""))
  if (any(init$w < 0)) stop("init$w may not be negative")  
  for (ff in 1:dd$nx_w) init$w[((ff - 1) * init$K + 1):(ff * init$K)] <- init$w[((ff - 1) * init$K + 1):(ff * init$K)] / sum(init$w[((ff - 1) * init$K + 1):(ff * init$K)])

  
  ##### init:  mu
  ##### ----------------------------------------------------  
  if (is.na(imu)){
    if (dd$p == 1){
      dist <- zR/(init$K + 1)
      init$mu <- seq(zMin+dist, zMax-dist, length=init$K)
    }else{
      dist <- zR/(init$K + 1)
      init$mu <- matrix(NA, nrow=init$K, ncol=dd$p)
      for (j in 1:dd$p) init$mu[,j] <- seq(zMin[j]+dist[j], zMax[j]-dist[j], length=init$K)
    }  
  }
  if (any(is.na(init$mu))) stop("NA in init$mu")          
  if (dd$p == 1){
    init$mu <- as.numeric(init$mu)
    if (length(init$mu) == CKmax & CKmax > init$K) init$mu <- init$mu[1:init$K]          
    if (length(init$mu) != init$K) stop(paste("init$mu must be of length ", init$K, sep=""))
    names(init$mu) <- paste("mu", 1:init$K, sep="")    
  }else{
    if (!is.matrix(init$mu)) stop("init$mu must be a matrix")
    if (ncol(init$mu) != dd$p) stop(paste("init$mu must have ", dd$p, " columns", sep=""))
    if (nrow(init$mu) != init$K) stop(paste("init$mu must have ", init$K, " rows", sep=""))
    rownames(init$mu) <- paste("j", 1:init$K, sep="")
    colnames(init$mu) <- paste("m", 1:dd$p, sep="")        
  }
  
  ##### init:  Sigma and Li
  ##### ----------------------------------------------------
  if (is.na(iSigma)){    
    if (is.na(iLi)){       ### Sigma and Li are computed from the data
      if (dd$p == 1){
        init$Sigma <- rep(zVar, init$K)
        names(init$Sigma) <- paste("Sigma", 1:init$K, sep="")
        init$Li <- sqrt(1 / init$Sigma)
        names(init$Li) <- paste("Li", 1:init$K, sep="")
      }else{
        init$Sigma <- matrix(rep(t(zVar), init$K), ncol=dd$p, byrow=TRUE)       ### initial Sigma matrix put in a big (K*p) x p matrix
        Sigmainv <- chol(zVar)        
        Sigmainv <- chol2inv(Sigmainv)
        Litmp <- t(chol(Sigmainv))
        init$Li <- rep(Litmp[lower.tri(Litmp, diag=TRUE)], init$K)
        rownames(init$Sigma) <- paste("j", rep(1:init$K, each=dd$p), ".", rep(1:dd$p, init$K), sep="")
        colnames(init$Sigma) <- paste("m", 1:dd$p, sep="")                
        names(init$Li) <- paste("Li", rep(1:init$K, each=dd$LTp), rep(dd$naamLTp, init$K), sep="")
      }        
    }else{                 ### Li is checked and Sigma is computed from Li
      if (any(is.na(init$Li))) stop("NA in init$Li")                    
      if (dd$p == 1){
        if (length(init$Li) == 1) init$Li <- rep(init$Li, init$K)
        if (length(init$Li) == CKmax & CKmax > init$K) init$Li <- init$Li[1:init$K]
        if (length(init$Sigma) != init$K) stop(paste("init$Sigma must be of length ", init$K, sep=""))
        init$Li <- as.numeric(init$Li)
        names(init$Li) <- paste("Li", 1:init$K, sep="")      
        if (any(init$Li <= 0)) stop("init$Li must be positive")
        init$Sigma <- (1 / init$Li)^2
        names(init$Sigma) <- paste("Sigma", 1:init$K, sep="")      
      }else{
        if (length(init$Li) == dd$LTp){
          tmpSigma <- matrix(0, nrow=dd$p, ncol=dd$p)
          tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init$Li
          tmpSigma <- tmpSigma %*% t(tmpSigma)
          err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
          if (class(err) == "try-error") stop("init$Li does not lead to a positive definite matrix")
          tmpSigma <- chol2inv(tmpSigma)
          init$Sigma <- matrix(rep(t(tmpSigma), init$K), ncol=dd$p, byrow=TRUE)
          init$Li <- rep(init$Li, init$K)
        }else{
          if (length(init$Li) == CKmax*dd$LTp & CKmax > init$K) init$Li <- init$Li[1:(init$K*dd$LTp)]
          if (length(init$Li) != init$K*dd$LTp) stop(paste("init$Li must be of length ", init$K*dd$LTp, sep=""))
          init$Sigma <- matrix(NA, ncol=dd$p, nrow=dd$p*init$K)
          for (j in 1:init$K){
            tmpSigma <- matrix(0, nrow=dd$p, ncol=dd$p)
            tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init$Li[((j-1)*dd$LTp+1):(j*dd$LTp)]
            tmpSigma <- tmpSigma %*% t(tmpSigma)
            err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
            if (class(err) == "try-error") stop(paste("the ", j,"-th block of init$Li does not lead to a positive definite matrix", sep=""))
            tmpSigma <- chol2inv(tmpSigma)
            init$Sigma[((j-1)*dd$p):(j*dd$p),] <- tmpSigma
          }  
        }
        rownames(init$Sigma) <- paste("j", rep(1:init$K, each=dd$p), ".", rep(1:dd$p, init$K), sep="")
        colnames(init$Sigma) <- paste("m", 1:dd$p, sep="")        
        names(init$Li) <- paste("Li", rep(1:init$K, each=dd$LTp), rep(dd$naamLTp, init$K), sep="")
      }  
    }  
  }else{                   ### Sigma is checked and Li is computed from Sigma
    if (any(is.na(init$Sigma))) stop("NA in init$Sigma")              
    if (dd$p == 1){
      if (length(init$Sigma) == 1) init$Sigma <- rep(init$Sigma, init$K)
      if (length(init$Sigma) == CKmax & CKmax > init$K) init$Sigma <- init$Sigma[1:init$K]      
      if (length(init$Sigma) != init$K) stop(paste("init$Sigma must be of length ", init$K, sep=""))
      init$Sigma <- as.numeric(init$Sigma)
      names(init$Sigma) <- paste("Sigma", 1:init$K, sep="")      
      if (any(init$Sigma <= 0)) stop("init$Sigma must be positive")
      init$Li <- sqrt(1 / init$Sigma)
      names(init$Li) <- paste("Li", 1:init$K, sep="")      
    }else{
      if (!is.matrix(init$Sigma)) stop("init$Sigma must be a matrix")
      if (ncol(init$Sigma) != dd$p) stop(paste("init$Sigma must have ", dd$p, " columns", sep=""))
      if (nrow(init$Sigma) == dd$p){
        if (any(init$Sigma[lower.tri(init$Sigma)] != t(init$Sigma)[lower.tri(init$Sigma)])) stop("init$Sigma must be a symmetric matrix")
        err <- try(Sigmainv <- chol(init$Sigma), silent=TRUE)
        if (class(err) == "try-error") stop("Cholesky decomposition of init$Sigma failed")
        Sigmainv <- chol2inv(Sigmainv)
        Litmp <- t(chol(Sigmainv))
        init$Li <- rep(Litmp[lower.tri(Litmp, diag=TRUE)], init$K)
      }else{
        if (nrow(init$Sigma) == CKmax*dd$p & CKmax > init$K) init$Sigma <- init$Sigma[1:(init$K*dd$p),]
        if (nrow(init$Sigma) != init$K*dd$p) stop(paste("init$Sigma must have ", init$K, " times ", dd$p, " rows", sep=""))
        init$Li <- numeric(0)
        for (j in 1:init$K){
          Sigmainv <- init$Sigma[((j-1)*dd$p+1):(j*dd$p),]
          if (any(Sigmainv[lower.tri(Sigmainv)] != t(Sigmainv)[lower.tri(Sigmainv)])) stop(paste(j, "-th block of init$Sigma is not symmetric", sep=""))
          err <- try(Sigmainv <- chol(Sigmainv), silent=TRUE)
          if (class(err) == "try-error") stop(paste("Cholesky decomposition of the ", j, "-th block of init$Sigma failed", sep=""))
          Sigmainv <- chol2inv(Sigmainv)
          Litmp <- t(chol(Sigmainv))
          init$Li <- c(init$Li, Litmp[lower.tri(Litmp, diag=TRUE)])
        }
      }
      rownames(init$Sigma) <- paste("j", rep(1:init$K, each=dd$p), ".", rep(1:dd$p, init$K), sep="")
      colnames(init$Sigma) <- paste("m", 1:dd$p, sep="")              
      names(init$Li) <- paste("Li", rep(1:init$K, each=dd$LTp), rep(dd$naamLTp, init$K), sep="")
    }  
  }    

  ##### init:  add Q to it
  ##### ----------------------------------------------------  
  if (dd$p == 1){
    init$Q <- 1 / init$Sigma
    names(init$Q) <- paste("Q", 1:init$K, sep="")    
  }else{
    init$Q <- numeric(0)
    for (j in 1:init$K){
      tmpLi <- diag(dd$p)
      tmpLi[lower.tri(tmpLi, diag=TRUE)] <- init$Li[((j-1)*dd$LTp + 1):(j*dd$LTp)]
      tmpQ <- tmpLi %*% t(tmpLi)
      init$Q <- c(init$Q, tmpQ[lower.tri(tmpQ, diag=TRUE)]) 
    }
    names(init$Q) <- paste("Q", rep(1:init$K, each=dd$LTp), rep(dd$naamLTp, init$K), sep="")
  }  
  
  ##### init:  gammaInv
  ##### ----------------------------------------------------  
  if (is.na(igammaInv)){
    if (dd$p == 1) init$gammaInv <- Czeta * zVar
    else           init$gammaInv <- Czeta * diag(zVar)
  }
  init$gammaInv <- as.numeric(init$gammaInv)
  if (length(init$gammaInv) == 1) init$gammaInv <- rep(init$gammaInv, dd$p)
  if (length(init$gammaInv) != dd$p) stop(paste("init$gammaInv must be of length ", dd$p, sep=""))
  if (any(is.na(init$gammaInv))) stop("NA in init$gammaInv")
  names(init$gammaInv) <- paste("gammaInv", 1:dd$p, sep="")

  
  ##### init:  r
  ##### ----------------------------------------------------
  if (is.na(ir)) init$r <- NMixMCMCinitr(z=initz, K=init$K, w=init$w[1:init$K], mu=init$mu, Sigma=init$Sigma, p=dd$p, n=dd$n)
  else           init$r <- NMixMCMCinitr(z=initz, K=init$K, w=init$w[1:init$K], mu=init$mu, Sigma=init$Sigma, p=dd$p, n=dd$n, initr=init$r) 
  rm(list="initz")
  

  ########## ========== Initials for chain 2 (if needed) ========== ##########
  ########## ====================================================== ##########
  if (PED){

    if (missing(init2)) init2 <- list()
    if (!is.list(init2)) stop("init2 must be a list")

    ininit2    <- names(init2)
    iy2        <- match("y",        ininit2, nomatch=NA)    
    iK2        <- match("K",        ininit2, nomatch=NA)
    iw2        <- match("w",        ininit2, nomatch=NA)
    imu2       <- match("mu",       ininit2, nomatch=NA)
    iSigma2    <- match("Sigma",    ininit2, nomatch=NA)
    iLi2       <- match("Li",       ininit2, nomatch=NA)  
    igammaInv2 <- match("gammaInv", ininit2, nomatch=NA)
    ir2        <- match("r",        ininit2, nomatch=NA)        
    
    ##### init2:  y
    ##### ----------------------------------------------------
    if (is.na(iy2)){
      init2$y <- NMixMCMCinity(y0=dd$y0, y1=dd$y1, censor=dd$censor, sd.init=sd(tmpinity),
                               are.Censored=dd$are.Censored, are.Right=dd$are.Right, are.Exact=dd$are.Exact, are.Left=dd$are.Left, are.Interval=dd$are.Interval,
                               p=dd$p, n=dd$n, random=TRUE)    
    }else{
      init2$y <- NMixMCMCinity(y0=dd$y0, y1=dd$y1, censor=dd$censor, sd.init=sd(tmpinity),
                               are.Censored=dd$are.Censored, are.Right=dd$are.Right, are.Exact=dd$are.Exact, are.Left=dd$are.Left, are.Interval=dd$are.Interval,
                               p=dd$p, n=dd$n, inity=init2$y)    
    }  


    ##### init2:  z (shifted and scaled y)
    ##### ----------------------------------------------------    
    initz2 <- (init2$y - matrix(rep(scale$shift, dd$n), ncol=dd$p, byrow=TRUE))/matrix(rep(scale$scale, dd$n), ncol=dd$p, byrow=TRUE)

    if (dd$p == 1){
      initz2 <- (init2$y - scale$shift)/scale$scale
      zBar2 <- mean(initz2)
      zMin2 <- min(initz2)
      zMax2 <- max(initz2)
    }else{
      initz2 <- (init2$y - matrix(rep(scale$shift, dd$n), ncol=dd$p, byrow=TRUE))/matrix(rep(scale$scale, dd$n), ncol=dd$p, byrow=TRUE)
      zBar2 <- apply(initz2, 2, mean)    ### will be used to determine prior$xi if this not given by the user
      zMin2 <- apply(initz2, 2, min)     
      zMax2 <- apply(initz2, 2, max)
    }
    zBar2[abs(zBar2) < 1e-14] <- 0  
    zVar2 <- var(initz2)                 ### will be used to determine init$Sigma and init$gammaInv if these not given by the user
    zVar2[abs(zVar2 - 1) < 1e-14] <- 1  
    zR2 <- zMax2 - zMin2                  ### will be used to determine init2$mu if this not given by the user
    zMid2 <- 0.5*(zMin2 + zMax2)    
    
    
    ##### init2:  K
    ##### ----------------------------------------------------  
    if (is.na(iK2)){
      if (prior$priorK == "fixed") init2$K <- CKmax
      else                         init2$K <- min(c(2, CKmax))
    }  
    if (length(init2$K) != 1) stop("init2$K must be of length 1")
    if (is.na(init2$K)) stop("NA in init2$K")
    if (init2$K <= 0 | init2$K > CKmax) stop("init2$K out of the range")

  
    ##### init2:  w
    ##### ----------------------------------------------------  
    if (is.na(iw2)){
      init2$w <- rDirichlet(1, rep(Cdelta, init2$K))
      if (dd$nx_w > 1) for (ff in 2:dd$nx_w) init2$w <- c(init2$w, rDirichlet(1, rep(Cdelta, init2$K)))
    }  
    init2$w <- as.numeric(init2$w)
    if (length(init2$w) == CKmax * dd$nx_w & CKmax > init2$K) init2$w <- init2$w[1:(init2$K * dd$nx_w)]
    if (dd$nx_w > 1) names(init2$w) <- paste("w", rep(1:init2$K, dd$nx_w), ".", rep(1:dd$nx_w, each = init2$K), sep="") else names(init2$w) <- paste("w", 1:init2$K, sep="")
    if (any(is.na(init2$w))) stop("NA in init2$w")  
    if (length(init2$w) != init2$K * dd$nx_w) stop(paste("init2$w must be of length ", init2$K * dd$nx_w, sep=""))
    if (any(init2$w < 0)) stop("init2$w may not be negative")  
    for (ff in 1:dd$nx_w) init2$w[((ff - 1) * init2$K + 1):(ff * init2$K)] <- init2$w[((ff - 1) * init2$K + 1):(ff * init2$K)] / sum(init2$w[((ff - 1) * init2$K + 1):(ff * init2$K)])

    
    ##### init2:  mu
    ##### ----------------------------------------------------  
    if (is.na(imu2)){
      tmpsd <- apply(tmpinitz, 2, sd)/init2$K    ## vector of length p      
      if (dd$p == 1){
        dist <- (zMax2 - zMin2)/(init2$K + 1)
        tmpxi <- seq(zMin2+dist, zMax2-dist, length=init2$K)
        init2$mu <- rnorm(init2$K, mean=tmpxi, sd=tmpsd)
        init2$mu <- init2$mu[order(init2$mu)]
      }else{
        dist <- (zMax2 - zMin2)/(init2$K + 1)
        init2$mu <- matrix(NA, nrow=init2$K, ncol=dd$p)
        for (j in 1:dd$p){
          tmpxi <- seq(zMin2[j]+dist[j], zMax2[j]-dist[j], length=init2$K)
          init2$mu[,j] <- rnorm(init2$K, mean=tmpxi, sd=tmpsd[j])
          init2$mu[,j] <- init2$mu[,j][order(init2$mu[,j])]
        }
      }  
    }    
    if (any(is.na(init2$mu))) stop("NA in init2$mu")          
    if (dd$p == 1){
      init2$mu <- as.numeric(init2$mu)
      if (length(init2$mu) == CKmax & CKmax > init2$K) init2$mu <- init2$mu[1:init2$K]          
      if (length(init2$mu) != init2$K) stop(paste("init2$mu must be of length ", init2$K, sep=""))
      names(init2$mu) <- paste("mu", 1:init2$K, sep="")    
    }else{
      if (!is.matrix(init2$mu)) stop("init2$mu must be a matrix")
      if (ncol(init2$mu) != dd$p) stop(paste("init2$mu must have ", dd$p, " columns", sep=""))
      if (nrow(init2$mu) != init2$K) stop(paste("init2$mu must have ", init2$K, " rows", sep=""))
      rownames(init2$mu) <- paste("j", 1:init2$K, sep="")
      colnames(init2$mu) <- paste("m", 1:dd$p, sep="")        
    }

    
    ##### init2:  Sigma and Li
    ##### ----------------------------------------------------  
    if (is.na(iSigma2)){    
      if (is.na(iLi2)){       ### Sigma and Li are computed from the data
        ctmp <- runif(init2$K, 0.1, 1.1)        
        if (dd$p == 1){
          init2$Sigma <- ctmp*zVar2
          names(init2$Sigma) <- paste("Sigma", 1:init2$K, sep="")
          init2$Li <- sqrt(1 / init2$Sigma)
          names(init2$Li) <- paste("Li", 1:init2$K, sep="")      
        }else{
          init2$Sigma <- ctmp[1]*zVar2
          Sigmainv <- chol(init2$Sigma)        
          Sigmainv <- chol2inv(Sigmainv)
          Litmp <- t(chol(Sigmainv))
          init2$Li <- Litmp[lower.tri(Litmp, diag=TRUE)]                    
          if (init2$K > 1){
            for (k in 2:init2$K){
              Sigmatmp <- ctmp[k]*zVar2
              init2$Sigma <- rbind(init2$Sigma, Sigmatmp)
              Sigmainv <- chol(Sigmatmp)        
              Sigmainv <- chol2inv(Sigmainv)
              Litmp <- t(chol(Sigmainv))
              init2$Li <- c(init2$Li, Litmp[lower.tri(Litmp, diag=TRUE)])              
            }
          }  
          rownames(init2$Sigma) <- paste("j", rep(1:init2$K, each=dd$p), ".", rep(1:dd$p, init2$K), sep="")
          colnames(init2$Sigma) <- paste("m", 1:dd$p, sep="")                
          names(init2$Li) <- paste("Li", rep(1:init2$K, each=dd$LTp), rep(dd$naamLTp, init2$K), sep="")
        }        
      }else{                 ### Li is checked and Sigma is computed from Li
        if (any(is.na(init2$Li))) stop("NA in init2$Li")                    
        if (dd$p == 1){
          if (length(init2$Li) == 1) init2$Li <- rep(init2$Li, init2$K)
          if (length(init2$Li) == CKmax & CKmax > init2$K) init2$Li <- init2$Li[1:init2$K]
          if (length(init2$Sigma) != init2$K) stop(paste("init2$Sigma must be of length ", init2$K, sep=""))
          init2$Li <- as.numeric(init2$Li)
          names(init2$Li) <- paste("Li", 1:init2$K, sep="")      
          if (any(init2$Li <= 0)) stop("init2$Li must be positive")
          init2$Sigma <- (1 / init2$Li)^2
          names(init2$Sigma) <- paste("Sigma", 1:init2$K, sep="")      
        }else{
          if (length(init2$Li) == dd$LTp){
            tmpSigma <- matrix(0, nrow=dd$p, ncol=dd$p)
            tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init2$Li
            tmpSigma <- tmpSigma %*% t(tmpSigma)
            err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
            if (class(err) == "try-error") stop("init2$Li does not lead to a positive definite matrix")
            tmpSigma <- chol2inv(tmpSigma)
            init2$Sigma <- matrix(rep(t(tmpSigma), init2$K), ncol=dd$p, byrow=TRUE)
            init2$Li <- rep(init2$Li, init2$K)
          }else{
            if (length(init2$Li) == CKmax*dd$LTp & CKmax > init2$K) init2$Li <- init2$Li[1:(init2$K*dd$LTp)]
            if (length(init2$Li) != init2$K*dd$LTp) stop(paste("init2$Li must be of length ", init2$K*dd$LTp, sep=""))
            init2$Sigma <- matrix(NA, ncol=dd$p, nrow=dd$p*init2$K)
            for (j in 1:init2$K){
              tmpSigma <- matrix(0, nrow=dd$p, ncol=dd$p)
              tmpSigma[lower.tri(tmpSigma, diag=TRUE)] <- init2$Li[((j-1)*dd$LTp+1):(j*dd$LTp)]
              tmpSigma <- tmpSigma %*% t(tmpSigma)
              err <- try(tmpSigma <- chol(tmpSigma), silent=TRUE)
              if (class(err) == "try-error") stop(paste("the ", j,"-th block of init2$Li does not lead to a positive definite matrix", sep=""))
              tmpSigma <- chol2inv(tmpSigma)
              init2$Sigma[((j-1)*dd$p):(j*dd$p),] <- tmpSigma
            }  
          }
          rownames(init2$Sigma) <- paste("j", rep(1:init2$K, each=dd$p), ".", rep(1:dd$p, init2$K), sep="")
          colnames(init2$Sigma) <- paste("m", 1:dd$p, sep="")        
          names(init2$Li) <- paste("Li", rep(1:init2$K, each=dd$LTp), rep(dd$naamLTp, init2$K), sep="")
        }  
      }  
    }else{                   ### Sigma is checked and Li is computed from Sigma
      if (any(is.na(init2$Sigma))) stop("NA in init2$Sigma")              
      if (dd$p == 1){
        if (length(init2$Sigma) == 1) init2$Sigma <- rep(init2$Sigma, init2$K)
        if (length(init2$Sigma) == CKmax & CKmax > init2$K) init2$Sigma <- init2$Sigma[1:init2$K]      
        if (length(init2$Sigma) != init2$K) stop(paste("init2$Sigma must be of length ", init2$K, sep=""))
        init2$Sigma <- as.numeric(init2$Sigma)
        names(init2$Sigma) <- paste("Sigma", 1:init2$K, sep="")      
        if (any(init2$Sigma <= 0)) stop("init2$Sigma must be positive")
        init2$Li <- sqrt(1 / init2$Sigma)
        names(init2$Li) <- paste("Li", 1:init2$K, sep="")      
      }else{
        if (!is.matrix(init2$Sigma)) stop("init2$Sigma must be a matrix")
        if (ncol(init2$Sigma) != dd$p) stop(paste("init2$Sigma must have ", dd$p, " columns", sep=""))
        if (nrow(init2$Sigma) == dd$p){
          if (any(init2$Sigma[lower.tri(init2$Sigma)] != t(init2$Sigma)[lower.tri(init2$Sigma)])) stop("init2$Sigma must be a symmetric matrix")
          err <- try(Sigmainv <- chol(init2$Sigma), silent=TRUE)
          if (class(err) == "try-error") stop("Cholesky decomposition of init2$Sigma failed")
          Sigmainv <- chol2inv(Sigmainv)
          Litmp <- t(chol(Sigmainv))
          init2$Li <- rep(Litmp[lower.tri(Litmp, diag=TRUE)], init2$K)
        }else{
          if (nrow(init2$Sigma) == CKmax*dd$p & CKmax > init2$K) init2$Sigma <- init2$Sigma[1:(init2$K*dd$p),]
          if (nrow(init2$Sigma) != init2$K*dd$p) stop(paste("init2$Sigma must have ", init2$K, " times ", dd$p, " rows", sep=""))
          init2$Li <- numeric(0)
          for (j in 1:init2$K){
            Sigmainv <- init2$Sigma[((j-1)*dd$p+1):(j*dd$p),]
            if (any(Sigmainv[lower.tri(Sigmainv)] != t(Sigmainv)[lower.tri(Sigmainv)])) stop(paste(j, "-th block of init2$Sigma is not symmetric", sep=""))
            err <- try(Sigmainv <- chol(Sigmainv), silent=TRUE)
            if (class(err) == "try-error") stop(paste("Cholesky decomposition of the ", j, "-th block of init2$Sigma failed", sep=""))
            Sigmainv <- chol2inv(Sigmainv)
            Litmp <- t(chol(Sigmainv))
            init2$Li <- c(init2$Li, Litmp[lower.tri(Litmp, diag=TRUE)])
          }
        }
        rownames(init2$Sigma) <- paste("j", rep(1:init2$K, each=dd$p), ".", rep(1:dd$p, init2$K), sep="")
        colnames(init2$Sigma) <- paste("m", 1:dd$p, sep="")              
        names(init2$Li) <- paste("Li", rep(1:init2$K, each=dd$LTp), rep(dd$naamLTp, init2$K), sep="")
      }  
    }    
  
    ##### init2:  gammaInv
    ##### ----------------------------------------------------  
    if (is.na(igammaInv2)){
      if (dd$p == 1) init2$gammaInv <- Czeta * zVar2 * runif(1, 0, dd$p)
      else           init2$gammaInv <- Czeta * diag(zVar2) * runif(dd$p, 0, dd$p)
    }
    init2$gammaInv <- as.numeric(init2$gammaInv)
    if (length(init2$gammaInv) == 1) init2$gammaInv <- rep(init2$gammaInv, dd$p)
    if (length(init2$gammaInv) != dd$p) stop(paste("init2$gammaInv must be of length ", dd$p, sep=""))
    if (any(is.na(init2$gammaInv))) stop("NA in init2$gammaInv")
    names(init2$gammaInv) <- paste("gammaInv", 1:dd$p, sep="")

  
    ##### init2:  r
    ##### ----------------------------------------------------
    if (is.na(ir2)) init2$r <- NMixMCMCinitr(z=initz2, K=init2$K, w=init2$w[1:init2$K], mu=init2$mu, Sigma=init2$Sigma, p=dd$p, n=dd$n)
    else            init2$r <- NMixMCMCinitr(z=initz2, K=init2$K, w=init2$w[1:init2$K], mu=init2$mu, Sigma=init2$Sigma, p=dd$p, n=dd$n, initr=init2$r) 
    rm(list="initz2")      
  }else{  ## end of if (PED)
    init2 <- NULL
  }  

  
  ########## ========== Values for the reversible jump MCMC ========== ##########
  ########## ========================================================= ##########
  if (missing(RJMCMC)) RJMCMC <- list()
  if (!is.list(RJMCMC)) stop("RJMCMC must be a list")

  inRJMCMC <- names(RJMCMC)
  iPaction <- match("Paction", ininit, nomatch=NA)
  iPsplit  <- match("Psplit",  ininit, nomatch=NA)
  iPbirth  <- match("Pbirth",  ininit, nomatch=NA)
  ipar.u1  <- match("par.u1",  ininit, nomatch=NA)
  ipar.u2  <- match("par.u2",  ininit, nomatch=NA)
  ipar.u3  <- match("par.u3",  ininit, nomatch=NA)  

  ##### RJMCMC:  Paction
  ##### ----------------------------------------------------
  if (is.na(iPaction)){
    actionAll <- 1
    RJMCMC$Paction <- c(1, 1, 1)/3
    #  RJMCMC$Paction <- c(0.2, 0.7, 0.1)     ### Dellaportas and Papageorgiou (2006)
  }   
  else{
    if (is.null(RJMCMC$Paction)){
      actionAll <- 1
      RJMCMC$Paction <- c(1, 1, 1)/3
      #  RJMCMC$Paction <- c(0.2, 0.7, 0.1)     ### Dellaportas and Papageorgiou (2006)      
    }else{
      actionAll <- 0
    }  
  }  
  if (length(RJMCMC$Paction) != 3) stop("RJMCMC$Paction must be of length 3")
  if (any(RJMCMC$Paction < 0)) stop("RJMCMC$Paction must be all non-negative")
  RJMCMC$Paction <- RJMCMC$Paction / sum(RJMCMC$Paction)
  CPaction <- as.numeric(RJMCMC$Paction)
  names(CPaction) <- names(RJMCMC$Paction) <- c("P.Gibbs.K", "P.split.combine", "P.birth.death")  
    
  ##### RJMCMC:  Psplit
  ##### ----------------------------------------------------  
  if (is.na(iPsplit)){
    if (CKmax == 1) RJMCMC$Psplit <- 0
    else            RJMCMC$Psplit <- c(1, rep(0.5, CKmax - 2), 0)
  }
  if (length(RJMCMC$Psplit) != CKmax) stop(paste("RJMCMC$Psplit must be of length ", CKmax, sep=""))
  if (any(RJMCMC$Psplit < 0)) stop("RJMCMC$Psplit must be all non-negative")
  if (any(RJMCMC$Psplit > 1)) stop("RJMCMC$Psplit must be all at most 1")
  if (RJMCMC$Psplit[CKmax] != 0) stop(paste("RJMCMC$Psplit[", CKmax, "] must be zero"))  
  if (CKmax > 1){
    if (RJMCMC$Psplit[1] != 1) stop(paste("RJMCMC$Psplit[", 1, "] must be one"))  
  }
  CPsplit <- as.numeric(RJMCMC$Psplit)
  names(CPsplit) <- names(RJMCMC$Psplit) <- paste("Psplit.", 1:CKmax, sep="")
    
  ##### RJMCMC:  Pbirth
  ##### ----------------------------------------------------  
  if (is.na(iPbirth)){
    if (CKmax == 1) RJMCMC$Pbirth <- 0
    else            RJMCMC$Pbirth <- c(1, rep(0.5, CKmax - 2), 0)
  }
  if (length(RJMCMC$Pbirth) != CKmax) stop(paste("RJMCMC$Pbirth must be of length ", CKmax, sep=""))
  if (any(RJMCMC$Pbirth < 0)) stop("RJMCMC$Pbirth must be all non-negative")
  if (any(RJMCMC$Pbirth > 1)) stop("RJMCMC$Pbirth must be all at most 1")
  if (RJMCMC$Pbirth[CKmax] != 0) stop(paste("RJMCMC$Pbirth[", CKmax, "] must be zero"))  
  if (CKmax > 1){
    if (RJMCMC$Pbirth[1] != 1) stop(paste("RJMCMC$Pbirth[", 1, "] must be one"))  
  }
  CPbirth <- as.numeric(RJMCMC$Pbirth)
  names(CPbirth) <- names(RJMCMC$Pbirth) <- paste("Pbirth.", 1:CKmax, sep="")
  
  ##### RJMCMC:  par.u1
  ##### ----------------------------------------------------  
  if (is.na(ipar.u1)) RJMCMC$par.u1 <- c(2, 2)
  if (length(RJMCMC$par.u1) != 2) stop("RJMCMC$par.u1 must be of length 2")
  if (any(RJMCMC$par.u1 <= 0)) stop("RJMCMC$par.u1 must be all positive")
  Cpar.u1 <- as.numeric(RJMCMC$par.u1)
  names(Cpar.u1) <- names(RJMCMC$par.u1) <- c("u1.1", "u1.2")
  
  ##### RJMCMC:  par.u2
  ##### ----------------------------------------------------
  if (dd$p == 1){
    if (is.na(ipar.u2)) RJMCMC$par.u2 <- c(1, 2*dd$p)
    if (length(RJMCMC$par.u2) != 2) stop("RJMCMC$par.u2 must be of length 2")
    if (any(RJMCMC$par.u2 <= 0)) stop("RJMCMC$par.u2 must be all positive")
    Cpar.u2 <- as.numeric(RJMCMC$par.u2)
    names(Cpar.u2) <- names(RJMCMC$par.u2) <- c("u2.1", "u2.2")    
  }else{
    if (is.na(ipar.u2)) RJMCMC$par.u2 <- rbind(matrix(c(rep(-1, dd$p-1), rep(1, dd$p-1)), ncol=2), c(1, 2*dd$p))
    if (!is.matrix(RJMCMC$par.u2))   stop("RJMCMC$par.u2 must be a matrix")
    if (nrow(RJMCMC$par.u2) != dd$p) stop(paste("RJMCMC$par.u2 must have ", dd$p, " rows", sep=""))
    if (ncol(RJMCMC$par.u2) != 2)    stop(paste("RJMCMC$par.u2 must have ", 2, " columns", sep=""))
    rownames(RJMCMC$par.u2) <- paste("u2.", 1:dd$p, sep="")
    colnames(RJMCMC$par.u2) <- paste(c(1, 2))
    if (any(RJMCMC$par.u2[dd$p,] <= 0)) stop(paste("RJMCMC$par.u2[", dd$p, ",] must be all positive", sep=""))    
    if (any(RJMCMC$par.u2[-dd$p,1] >= RJMCMC$par.u2[-dd$p,2])) stop(paste("The first column of RJMCMC$par.u2[-", dd$p, ",] must be strictly lower than the second column", sep=""))
    Cpar.u2 <- as.numeric(t(RJMCMC$par.u2))
    names(Cpar.u2) <- paste("u2.", rep(1:dd$p, each=2), ".", rep(1:2, dd$p), sep="")    
  }  
    
  ##### RJMCMC:  par.u3
  ##### ----------------------------------------------------  
  if (dd$p == 1){
    if (is.na(ipar.u3)) RJMCMC$par.u3 <- c(1, dd$p)
    if (length(RJMCMC$par.u3) != 2) stop("RJMCMC$par.u3 must be of length 2")
    if (any(RJMCMC$par.u3 <= 0))    stop("RJMCMC$par.u3 must be all positive")
    Cpar.u3 <- as.numeric(RJMCMC$par.u3)
    names(Cpar.u3) <- names(RJMCMC$par.u3) <- c("u3.1", "u3.2")    
  }else{
    if (is.na(ipar.u3)) RJMCMC$par.u3 <- rbind(matrix(c(rep(0, dd$p-1), rep(1, dd$p-1)), ncol=2), c(1, dd$p))
    if (!is.matrix(RJMCMC$par.u3)) stop("RJMCMC$par.u3 must be a matrix")
    if (nrow(RJMCMC$par.u3) != dd$p) stop(paste("RJMCMC$par.u3 must have ", dd$p, " rows", sep=""))
    if (ncol(RJMCMC$par.u3) != 2)    stop(paste("RJMCMC$par.u3 must have ", 2, " columns", sep=""))
    rownames(RJMCMC$par.u3) <- paste("u3.", 1:dd$p, sep="")
    colnames(RJMCMC$par.u3) <- paste(c(1, 2))
    if (any(RJMCMC$par.u3[dd$p,] <= 0)) stop(paste("RJMCMC$par.u3[", dd$p, ",] must be all positive", sep=""))    
    if (any(RJMCMC$par.u3[-dd$p,1] >= RJMCMC$par.u3[-dd$p,2])) stop(paste("The first column of RJMCMC$par.u3[-", dd$p, ",] must be strictly lower than the second column", sep=""))
    Cpar.u3 <- as.numeric(t(RJMCMC$par.u3))
    names(Cpar.u3) <- paste("u3.", rep(1:dd$p, each=2), ".", rep(1:2, dd$p), sep="")    
  }  

  ##### RJMCMC:  concetenate
  ##### ----------------------------------------------------
  CRJMCMC <- c(CPaction, CPsplit, CPbirth, Cpar.u1, Cpar.u2, Cpar.u3)

  ########## ========== RETURN if onlyInit ========== ##########
  ########## ======================================== ##########  
  if (onlyInit) return(list(prior=prior, init=init, init2=init2, scale=scale, RJMCMC=RJMCMC))
  

  ########## ========== nMCMC ========== ##########
  ########## =========================== ##########
  if (length(nMCMC) != 4) stop("nMCMC must be of length 4")
  if (is.null(names(nMCMC))) names(nMCMC) <- c("burn", "keep", "thin", "info")
  names.nMCMC <- names(nMCMC)
  if (!match("burn", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("burn"), "] must be specified", sep=""))
  else                                        n.burn <- nMCMC["burn"]
  if (!match("keep", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("keep"), "] must be specified", sep=""))
  else                                        n.keep <- nMCMC["keep"]
  if (!match("thin", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("thin"), "] must be specified", sep=""))
  else                                        n.thin <- nMCMC["thin"]
  if (!match("info", names.nMCMC, nomatch=0)) stop(paste("nMCMC[", dQuote("info"), "] must be specified", sep=""))
  else                                        n.info <- nMCMC["info"]
  nMCMC <- c(n.burn, n.keep, n.thin, n.info)
  names(nMCMC) <- c("burn", "keep", "thin", "info")  
  if (nMCMC["burn"] < 0) stop(paste("nMCMC[", dQuote("burn"), "] must be non-negative", sep=""))
  if (nMCMC["keep"] <= 0) stop(paste("nMCMC[", dQuote("keep"), "] must be positive", sep=""))  
  if (nMCMC["thin"] <= 0) stop(paste("nMCMC[", dQuote("thin"), "] must be positive", sep=""))
  if (nMCMC["info"] <= 0 | nMCMC["info"] > max(nMCMC["burn"], nMCMC["keep"])) nMCMC["info"] <- max(nMCMC["burn"], nMCMC["keep"])
  
#  RET <- "RESULT"
#  attr(RET, "prior") <- prior
#  attr(RET, "init") <- init
#  attr(RET, "init2") <- init2  
#  attr(RET, "RJMCMC") <- RJMCMC
#  attr(RET, "nMCMC") <- nMCMC
#
#  attr(RET, "Cprior") <- list(integer=Cinteger, double=Cdouble)
#  attr(RET, "CRJMCMC") <- CRJMCMC    


  ##### Parameters passed to C++ which will be stored also in the resulting object (to be able to use them in related functions)
  ##### ========================================================================================================================
  Cpar <- list(z0          = z0,
               z1          = z1,
               censor      = dd$censor,               
               dimy        = c(p=dd$p, n=dd$n),
               priorInt    = Cinteger,
               priorDouble = Cdouble,
               x_w         = as.integer(c(dd$nx_w, dd$x_w)))
  rm(list=c("Cinteger", "Cdouble"))
    
  
  ########## ========== Run MCMC ========== ##########
  ########## ============================== ##########
  if (PED){
    if (parallel){
      #require("parallel")

      if (parallel::detectCores() < 2) warning("It does not seem that at least 2 CPU cores are available needed for efficient parallel generation of the two chains.")
      if (missing(cltype)) cl <- parallel::makeCluster(2) else cl <- parallel::makeCluster(2, type = cltype)
      cat(paste("Parallel MCMC sampling of two chains started on ", date(), ".\n", sep=""))      
      RET <- parallel::parLapply(cl, 1:2, NMixMCMCwrapper,
                                 scale = scale, prior = prior, inits = list(init, init2), Cpar = Cpar, RJMCMC = RJMCMC, CRJMCMC = CRJMCMC,
                                 actionAll = actionAll, nMCMC = nMCMC, keep.chains = keep.chains, PED = TRUE,
                                 dens.zero = dens.zero, lx_w = dd$lx_w)
      cat(paste("Parallel MCMC sampling finished on ", date(), ".\n", sep=""))
      parallel::stopCluster(cl)      
    }else{
      RET <- lapply(1:2, NMixMCMCwrapper,
                         scale = scale, prior = prior, inits = list(init, init2), Cpar = Cpar, RJMCMC = RJMCMC, CRJMCMC = CRJMCMC,
                         actionAll = actionAll, nMCMC = nMCMC, keep.chains = keep.chains, PED = TRUE,
                         dens.zero = dens.zero, lx_w = dd$lx_w)
    }
    
    cat(paste("\nComputation of penalized expected deviance started on ", date(), ".\n", sep=""))
    
    if (prior$priorK == "fixed"){
      resPED <- .C("NMix_PED", PED               = double(5),
                               pm.indDevObs      = double(dd$n),
                               pm.indpopt        = double(dd$n),
                               pm.windpopt       = double(dd$n),
                               invalid.indDevObs = integer(dd$n),
                               invalid.indpopt   = integer(dd$n),
                               invalid.windpopt  = integer(dd$n),                   
                               sum.ISweight      = double(dd$n),
                               #ch.ISweight       = double(dd$n*nMCMC["keep"]),                   
                               err               = integer(1),
                               y0                = as.double(t(z0)),
                               y1                = as.double(t(z1)),
                               censor            = as.integer(t(dd$censor)),
                               nxw_xw            = as.integer(Cpar$x_w),                     
                               dimy              = as.integer(c(dd$p, dd$n)),
                               chK1              = as.integer(RET[[1]]$K),
                               chw1              = as.double(t(RET[[1]]$w)),
                               chmu1             = as.double(t(RET[[1]]$mu)),
                               chLi1             = as.double(t(RET[[1]]$Li)),
                               chK2              = as.integer(RET[[1]]$K),
                               chw2              = as.double(t(RET[[2]]$w)),
                               chmu2             = as.double(t(RET[[2]]$mu)),
                               chLi2             = as.double(t(RET[[2]]$Li)),
                               M                 = as.integer(nMCMC["keep"]),
                               Kmax              = as.integer(CKmax),
                               Krandom           = as.integer(0),
                               Dens.ZERO         = as.double(dens.zero),
                               EMin              = as.double(EMin),
                   PACKAGE = thispackage)
    }else{
      resPED <- .C("NMix_PED", PED               = double(5),
                               pm.indDevObs      = double(dd$n),
                               pm.indpopt        = double(dd$n),
                               pm.windpopt       = double(dd$n),
                               invalid.indDevObs = integer(dd$n),
                               invalid.indpopt   = integer(dd$n),
                               invalid.windpopt  = integer(dd$n),
                               sum.ISweight      = double(dd$n),
                               #ch.ISweight       = double(dd$n*nMCMC["keep"]),
                               err               = integer(1),
                               y0                = as.double(t(z0)),
                               y1                = as.double(t(z1)),
                               censor            = as.integer(t(dd$censor)),
                               nxw_xw            = as.integer(Cpar$x_w),                   
                               dimy              = as.integer(c(dd$p, dd$n)),
                               chK1              = as.integer(RET[[1]]$K),
                               chw1              = as.double(RET[[1]]$w),
                               chmu1             = as.double(RET[[1]]$mu),
                               chLi1             = as.double(RET[[1]]$Li),
                               chK2              = as.integer(RET[[2]]$K),
                               chw2              = as.double(RET[[2]]$w),
                               chmu2             = as.double(RET[[2]]$mu),
                               chLi2             = as.double(RET[[2]]$Li),
                               M                 = as.integer(nMCMC["keep"]),
                               Kmax              = as.integer(CKmax),
                               Krandom           = as.integer(1),
                               Dens.ZERO         = as.double(dens.zero),
                               EMin              = as.double(EMin),                   
                   PACKAGE = thispackage)
    }
    
    cat(paste("Computation of penalized expected deviance finished on ", date(), ".\n", sep=""))
    if (resPED$err) stop("Something went wrong.")

    names(resPED$PED) <- c("D.expect", "p(opt)", "PED", "wp(opt)", "wPED")
    RET$PED <- resPED$PED
    
    detS      <- prod(scale$scale)
    idetS     <- 1 / detS
    log.idetS <- -log(detS)

    RET$PED["D.expect"] <- RET$PED["D.expect"] - 2*dd$n*log.idetS     ## correction for scale transformation
    RET$PED["PED"]      <- RET$PED["PED"] - 2*dd$n*log.idetS          ## correction for scale transformation
    RET$PED["wPED"]     <- RET$PED["wPED"] - 2*dd$n*log.idetS         ## correction for scale transformation
    
    RET$D         <- resPED$pm.indDevObs - 2*log.idetS                        ## correction for scale transformation
    RET$popt      <- resPED$pm.indpopt
    RET$wpopt     <- resPED$pm.windpopt
    RET$inv.D     <- resPED$invalid.indDevObs
    RET$inv.popt  <- resPED$invalid.indpopt
    RET$inv.wpopt <- resPED$invalid.windpopt    
    RET$sumISw    <- resPED$sum.ISweight
    #RET$ISw       <- matrix(resPED$ch.ISweight, nrow=dd$n, byrow=TRUE)

    class(RET) <- "NMixMCMClist"
  }else{
    RET <- NMixMCMCwrapper(chain = 1,
                           scale = scale, prior = prior, inits = list(init), Cpar = Cpar, RJMCMC = RJMCMC, CRJMCMC = CRJMCMC,
                           actionAll = actionAll, nMCMC = nMCMC, keep.chains = keep.chains,
                           PED = FALSE, dens.zero = dens.zero, lx_w = dd$lx_w)
  }  
  
  return(RET)
}
