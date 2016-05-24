##########################################################################
## MCMCmnl.R samples from the posterior distribution of a multinomial
## logit model using Metropolis-Hastings.
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## KQ 12/22/2004
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

## parse formula and return a list that contains the model response
## matrix as element one, the model matrix as element two,
## the column names of X as element three, the rownames of
## X and y as element four, and the number of choices in the largest
## choice set in element five
"parse.formula.mnl" <- function(formula, data, baseline=NULL,
                                intercept=TRUE){
  
  ## extract Y, X, and variable names for model formula and frame
  mt <- terms(formula, data=data)
  if(missing(data)) data <- sys.frame(sys.parent())
  mf <- match.call(expand.dots = FALSE)
  mf$intercept <- mf$baseline <- NULL
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  
  mf <- eval(mf, sys.frame(sys.parent()))
  if (!intercept){
    attributes(mt)$intercept <- 0
  }
  
  
  ## deal with Y
  Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
  if (ncol(Y)==1){
    Y <- factor(Y)
    number.choices <- length(unique(Y))
    choice.names <- sort(unique(Y))
    Ymat <- matrix(NA, length(Y), number.choices)
    colnames(Ymat) <- choice.names
    for (i in 1:(number.choices)){
      Ymat[,i] <- as.numeric(Y==choice.names[i])          
    }
  }
  else{
    ## this block will allow for nonconstant choice sets
    number.choices <- ncol(Y)
    Ytemp <- Y
    Ytemp[Y== -999] <- NA
    if ( min(unique(array(Y)) %in% c(-999,0,1))==0 ||
        min(apply(Ytemp, 1, sum, na.rm=TRUE) == rep(1, nrow(Y)))==0){
      stop("Y is a matrix but it is not composed of 0/1/-999 values\n and/or rows do not sum to 1\n")
    }
    Ymat <- Y
    choice.names <- colnames(Y)
  }
  colnames(Ymat) <- choice.names
  rownames(Ymat) <- 1:nrow(Ymat)
  
  #Y.long <- matrix(t(Ymat), length(Ymat), 1)
  #colnames(Y.long) <- "Y"
  #rownames(Y.long) <- rep(1:nrow(Ymat), rep(length(choice.names), nrow(Ymat)))
  #rownames(Y.long) <- paste(rownames(Y.long), choice.names, sep=".")
  #group.id <- rep(1:nrow(Ymat), rep(ncol(Ymat), nrow(Ymat)))
  
  ## deal with X
  ## null model support
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
  X <- as.matrix(X)         # X matrix
  xvars <- dimnames(X)[[2]] # X variable names
  xobs  <- dimnames(X)[[1]] # X observation names
  
  if (is.null(baseline)){
    baseline <- choice.names[1]
  }
  if (! baseline %in% choice.names){
    stop("'baseline' not consistent with observed choice levels in y\n")
  }
  
  ## deal with choice specific covariates
  choicevar.indic <- rep(FALSE, length(xvars)) # indicators for choice
                                               # specific variables
  choicevar.indic[grep("^choicevar\\(", xvars)] <- TRUE
  if (sum(choicevar.indic) > 0){
    cvarname1.vec <- rep(NA, sum(choicevar.indic))
    cvarname2.vec <- rep(NA, sum(choicevar.indic))
    counter <- 0
    for (i in 1:length(xvars)){
      if (choicevar.indic[i]){
        counter <- counter + 1
        vn2 <- strsplit(xvars[i], ",")
        vn3 <- strsplit(vn2[[1]], "\\(")
        vn4 <- strsplit(vn3[[3]], "=")
        cvarname1 <- vn3[[2]][1]
        cvarname1 <- strsplit(cvarname1, "\"")[[1]]
        cvarname1 <- cvarname1[length(cvarname1)]
        cvarname2 <- vn4[[1]][length(vn4[[1]])]
        cvarname2 <- strsplit(cvarname2, "\"")[[1]][2]
        if (! cvarname2 %in% choice.names){
          stop("choicelevel that was set in choicevar() not consistent with\n observed choice levels in y")
        }
        cvarname1.vec[counter] <- cvarname1
        cvarname2.vec[counter] <- cvarname2
        xvars[i] <- paste(cvarname1, cvarname2, sep=".")
      }
    }
    
    X.cho <- X[, choicevar.indic]
    X.cho.mat <- matrix(NA, length(choice.names)*nrow(X),
                        length(unique(cvarname1.vec)))
    rownames(X.cho.mat) <- rep(rownames(X), rep(length(choice.names), nrow(X)))
    rownames(X.cho.mat) <- paste(rownames(X.cho.mat), choice.names, sep=".")
    colnames(X.cho.mat) <- unique(cvarname1.vec)
    choice.names.n <- rep(choice.names, nrow(X))
    for (j in 1:length(unique(cvarname1.vec))){
      for (i in 1:length(cvarname2.vec)){
        if (colnames(X.cho.mat)[j] == cvarname1.vec[i]){
          X.cho.mat[choice.names.n==cvarname2.vec[i], j] <-
            X.cho[,i]
        }
      }
    }
  }

  
  ## deal with individual specific covariates
  xvars.ind.mat <- rep(xvars[!choicevar.indic],
                       rep(length(choice.names),
                           sum(!choicevar.indic)))
  xvars.ind.mat <- paste(xvars.ind.mat, choice.names, sep=".")

  if (sum(!choicevar.indic) > 0){
    X.ind.mat <- X[,!choicevar.indic] %x% diag(length(choice.names))
    colnames(X.ind.mat) <- xvars.ind.mat
    rownames(X.ind.mat) <- rep(rownames(X), rep(length(choice.names), nrow(X)))
    rownames(X.ind.mat) <- paste(rownames(X.ind.mat), choice.names, sep=".")
    
    ## delete columns correpsonding to the baseline choice
    ivarname1 <- strsplit(xvars.ind.mat, "\\.")
    ivar.keep.indic <- rep(NA, ncol(X.ind.mat))
    for (i in 1:ncol(X.ind.mat)){
      ivar.keep.indic[i] <- ivarname1[[i]][length(ivarname1[[i]])] != baseline
    }
    X.ind.mat <- X.ind.mat[,ivar.keep.indic]
  }
  
  if (sum(choicevar.indic) > 0 & sum(!choicevar.indic) > 0){
    X <- cbind(X.cho.mat, X.ind.mat)
  }
  else if (sum(!choicevar.indic) > 0){
    X <- X.ind.mat
  }
  else if (sum(choicevar.indic) > 0){
    X <- X.cho.mat
  }
  else {
    stop("X matrix appears to have neither choice-specific nor individual-specific variables.\n")
  }
  #Y <- Y.long
  xvars <- colnames(X)
  xobs <- rownames(X)
  
  return(list(Ymat, X, xvars, xobs, number.choices))
  
}

## dummy function used to handle choice-specific covariates
"choicevar" <- function(var, varname, choicelevel){
  junk1 <- varname
  junk2 <- choicelevel
  return(var)
}


## MNL log-posterior function (used to get starting values)
## vector Y without NAs
"mnl.logpost.noNA" <- function(beta, new.Y, X, b0, B0){
  nobs <- length(new.Y)
  ncat <- nrow(X) / nobs
  Xb <- X %*% beta
  Xb <- matrix(Xb, byrow=TRUE, ncol=ncat)
  indices <- cbind(1:nobs, new.Y)
  Xb.reform <- Xb[indices]
  eXb <- exp(Xb)
  #denom <- log(apply(eXb, 1, sum))
  z <- rep(1, ncat)
  denom <- log(eXb %*% z)

  log.prior <- 0.5 * t(beta - b0) %*% B0 %*% (beta - b0)
  
  return(sum(Xb.reform - denom) + log.prior)
}

## MNL log-posterior function (used to get starting values)
## matrix Y with NAs
"mnl.logpost.NA" <- function(beta, Y, X, b0, B0){
    k <- ncol(X)
    numer <- exp(X %*% beta)
    numer[Y== -999] <- NA
    numer.mat <- matrix(numer, nrow(Y), ncol(Y), byrow=TRUE)
    denom <- apply(numer.mat, 1, sum, na.rm=TRUE)
    choice.probs <- numer.mat / denom
    Yna <- Y
    Yna[Y == -999] <- NA  
    log.like.mat <- log(choice.probs) * Yna
    log.like <- sum(apply(log.like.mat, 1, sum, na.rm=TRUE))

    log.prior <- 0.5 * t(beta - b0) %*% B0 %*% (beta - b0)

    return(log.like + log.prior)
}




"MCMCmnl" <-
  function(formula, baseline=NULL, data=NULL, 
           burnin = 1000, mcmc = 10000, thin=1,
           mcmc.method = c("IndMH", "RWM", "slice"),
           tune = 1.0, tdf=6, verbose = 0, seed = NA,
           beta.start = NA, b0 = 0, B0 = 0, ...) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)
    if (tdf <= 0){
      stop("degrees of freedom for multivariate-t proposal must be positive.\n Respecify tdf and try again.\n") 
    }
    
    ## seeds
    seeds <- form.seeds(seed) 
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]
    
    ## form response and model matrix
    holder <- parse.formula.mnl(formula=formula, baseline=baseline,
                                data=data)
    Y <- holder[[1]]
    ## check to make sure baseline category is always available in choiceset
    if (is.null(baseline)){
      if (max(Y[,1] == -999) == 1){
        stop("Baseline choice not available in all choicesets.\n Respecify baseline category and try again.\n")
      }
    }
    else{
      if (max(Y[,baseline] == -999) == 1){
        stop("Baseline choice not available in all choicesets.\n Respecify baseline category and try again.\n")
      }      
    }
    X <- holder[[2]]
    xnames <- holder[[3]]
    xobs <- holder[[4]]
    number.choices <- holder[[5]]
    K <- ncol(X)  # number of covariates
    
    ## form the tuning parameter
    tune <- vector.tune(tune, K)

    ## priors and starting values 
    mvn.prior <- form.mvn.prior(b0, B0, K)
    b0 <- mvn.prior[[1]]
    B0 <- mvn.prior[[2]]

    beta.init <- rep(0, K)
    cat("Calculating MLEs and large sample var-cov matrix.\n")
    cat("This may take a moment...\n")
    if (max(is.na(Y))){
      optim.out <- optim(beta.init, mnl.logpost.NA, method="BFGS",
                         control=list(fnscale=-1),
                         hessian=TRUE, Y=Y, X=X, b0=b0, B0=B0)
    }
    else{
      new.Y <- apply(Y==1, 1, which)
      optim.out <- optim(beta.init, mnl.logpost.noNA, method="BFGS",
                         control=list(fnscale=-1),
                         hessian=TRUE, new.Y=new.Y, X=X, b0=b0, B0=B0)      
    }
    cat("Inverting Hessian to get large sample var-cov matrix.\n")
    ##V <- solve(-1*optim.out$hessian)
    V <- chol2inv(chol(-1*optim.out$hessian))
    beta.mode <- matrix(optim.out$par, K, 1)
    
    
    if (is.na(beta.start) || is.null(beta.start)){
      beta.start <- matrix(optim.out$par, K, 1)
    }
    else if(is.null(dim(beta.start))) {
      beta.start <- matrix(beta.start, K, 1)
    }
    else if (length(beta.start != K)){
      stop("beta.start not of appropriate dimension\n")
    }
      
    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, dim(X)[2] )
    posterior <- NULL    

    if (mcmc.method=="RWM"){
      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCmnlMH",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       tune=tune, lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose),
                       betastart=beta.start, betamode=beta.mode,
                       b0=b0, B0=B0,
                       V=V, RW=as.integer(1), tdf=as.double(tdf)) 
      
      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCmnl Posterior Sample")
    }
    else if (mcmc.method=="IndMH"){
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCmnlMH",
                       sample.nonconst=sample, Y=Y, X=X,
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       tune=tune, lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose),
                       betastart=beta.start, betamode=beta.mode,
                       b0=b0, B0=B0,
                       V=V, RW=as.integer(0), tdf=as.double(tdf)) 

      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCmnl Posterior Sample")
      
    }
    else if (mcmc.method=="slice"){
      ## call C++ code to draw sample
      auto.Scythe.call(output.object="posterior", cc.fun.name="MCMCmnlslice",
                       sample.nonconst=sample, Y=Y, X=X, 
                       burnin=as.integer(burnin),
                       mcmc=as.integer(mcmc), thin=as.integer(thin),
                       lecuyer=as.integer(lecuyer),
                       seedarray=as.integer(seed.array),
                       lecuyerstream=as.integer(lecuyer.stream),
                       verbose=as.integer(verbose), betastart=beta.start,
                       b0=b0, B0=B0, V=V) 
      
      ## put together matrix and build MCMC object to return
      output <- form.mcmc.object(posterior, names=xnames,
                                 title="MCMCmnl Posterior Sample")

    }

    return(output) 
   
  }



