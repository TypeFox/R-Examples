.onAttach <-
  function (libname, pkgname) 
  {
    loadmsg <- "\nNote that there are changes to the names of some functions in version 1.0.3. See NEWS.\n"
    packageStartupMessage(loadmsg, domain = NULL, appendLF = TRUE)
  }

`randomLCA` <-
  function(patterns,freq=NULL,nclass=2,calcSE=TRUE,notrials=20,
           random=FALSE,byclass=FALSE,quadpoints=21,constload=TRUE,blocksize=dim(patterns)[2],
           level2=FALSE,probit=FALSE,level2size=blocksize,
           qniterations=5,penalty=0.001,verbose=FALSE,seed = as.integer(runif(1, 0, .Machine$integer.max))) {
    set.seed(seed)
    if (quadpoints > 190)
      stop("Maximum of 190 quadrature points\n")
    cl <- match.call()
    # check that patterns are either 0 or 1
    if (any(apply(as.matrix(patterns),2,function(x) any((x!=0)&(x!=1)&!is.na(x)))))
      stop("patterns must consist of either 0 or 1")    
    # check that patterns doesn't contain column which is all missing
    if (any(apply(as.matrix(patterns),2,function(x) all(is.na(x)))))
      stop("patterns cannot contain columns consisting entirely of missing")
    # check that freq are all >= 0
    if (!missing(freq)) {
      if (any(freq<0))
        stop("frequencies must be greater than or equal to zero")
    }
    if (random & ((dim(patterns)[2] %% blocksize)!=0))
      stop("number of outcomes must be a multiple of blocksize")
    # if no frequencies given, then assume that the data needs to be summarised
    if (missing(freq) | is.null(freq)) {
      pats <- apply(patterns, 1, function(x) {paste(ifelse(is.na(x),"N",x),collapse="")})
      tpats <- table(pats)
      freq <- as.numeric(tpats)
      newpatterns <- unlist(strsplit(names(tpats),split=""))
      newpatterns <- ifelse(newpatterns=="N",NA_character_,newpatterns)
      newpatterns <- as.data.frame(matrix(as.numeric(newpatterns),byrow=TRUE,ncol=dim(patterns)[2]))
      if (is.null(names(patterns))) names(newpatterns) <- paste("X",1:dim(patterns)[2],sep="")
      else names(newpatterns) <- names(patterns)
      patterns <- newpatterns
    }
    else {
      # check that freq doesn't contain missing
      if (any(is.na(freq))) stop("freq cannot contain missing values")
      # remove any observations with frequency of zero
#       patterns <- patterns[freq!=0,]
#       freq <- freq[freq!=0]
    }
# determine df
    nparams <- dim(patterns)[2]*nclass
    nparams <- nparams+nclass-1
    if (random) {
        if (level2) {
          nparams <- nparams+ifelse(constload,1,level2size)*ifelse(byclass,nclass,1)
          nparams <- nparams+ifelse(byclass,nclass,1)
        } else nparams <- nparams+ifelse(constload,1,min(blocksize,dim(patterns)[2]))*ifelse(byclass,nclass,1)
    }
    df <- 2^dim(patterns)[2]-nparams-1
    #print(paste('df = ',df))
    nonident <- FALSE
    if (df < 0) nonident <- TRUE
    if ((nclass==2) & (dim(patterns)[2]<3)) nonident <- TRUE
    if ((nclass==3) & (dim(patterns)[2]<5)) nonident <- TRUE
    if ((nclass==4) & (dim(patterns)[2]<5)) nonident <- TRUE
    if ((nclass==5) & (dim(patterns)[2]<5)) nonident <- TRUE
    if (nonident) stop("Model is not identifiable - decrease classes or random effects")
    if (!random) initmodel <- bestlca(patterns,freq=freq,nclass=nclass,
            calcSE=(calcSE & !random),notrials=notrials,probit=probit,penalty=penalty,verbose=verbose)
    else {
      if (!level2) {
        initmodel <- bestlca(patterns,freq=freq,nclass=nclass,
                             calcSE=FALSE,notrials=notrials,probit=probit,penalty=penalty,verbose=verbose)
        # work out how many lambda coefs there are
        if (constload) nlambda <- 1
        else nlambda <- min(dim(patterns)[2],blocksize)
        # now fit the simplest random efefcts model ie with constant loading
        initmodel <- fitAdaptRandom(patterns,freq=freq,
                                                nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
                                                initclassp=initmodel$classp,initlambdacoef=NULL,
                                                gh=norm.gauss.hermite(quadpoints),
                                                constload=TRUE,probit=probit,byclass=FALSE,qniterations=qniterations,
                                                penalty=penalty,verbose=verbose)
        # fit with variable loading if required
        if (!constload)  initmodel <- fitAdaptRandom(patterns,freq=freq,
                                                nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
                                                initclassp=initmodel$classp,
                                                 initlambdacoef=rep(initmodel$lambdacoef,nlambda),
                                                 gh=norm.gauss.hermite(quadpoints),
                                                  constload=constload,blocksize=blocksize,
                                                 probit=probit,byclass=FALSE,qniterations=qniterations,
                                                  penalty=penalty,verbose=verbose)
        if (byclass)  {
          initlambdacoef <- matrix(rep(initmodel$lambdacoef,nclass),nrow=nclass,byrow=TRUE)
          initmodel <- fitAdaptRandom(patterns,freq=freq,
                                                  nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
                                                  initclassp=initmodel$classp,
                                                  initlambdacoef=initlambdacoef,
                                                  gh=norm.gauss.hermite(quadpoints),
                                                  constload=constload,blocksize=blocksize,
                                                  probit=probit,byclass=byclass,qniterations=qniterations,
                                                  penalty=penalty,verbose=verbose)
        }
      } else {
        # 2 level models
        initmodel <- bestlca(patterns,freq=freq,nclass=nclass,
                             calcSE=FALSE,notrials=notrials,probit=probit,penalty=penalty,verbose=verbose)
        # work out how many lambda coefs there are
        if (constload) nlambda <- 1
        else nlambda <- min(dim(patterns)[2],blocksize)
        # now fit the simplest random efefcts model
        initmodel <- fitAdaptRandom(patterns,freq=freq,
                                                nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
                                                initclassp=initmodel$classp,initlambdacoef=NULL,
                                                gh=norm.gauss.hermite(quadpoints),
                                                constload=TRUE,probit=probit,byclass=FALSE,qniterations=qniterations,
                                                penalty=penalty,verbose=verbose)
        # fit with variable loading if required
        if (!constload)  initmodel <- fitAdaptRandom(patterns,freq=freq,
                                                                 nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
                                                                 initclassp=initmodel$classp,
                                                                 initlambdacoef=rep(initmodel$lambdacoef,nlambda),
                                                                 gh=norm.gauss.hermite(quadpoints),
                                                                 constload=constload,blocksize=level2size,
                                                                 probit=probit,byclass=FALSE,qniterations=qniterations,
                                                                 penalty=penalty,verbose=verbose)
        if (byclass)  {
          initlambdacoef <- matrix(rep(initmodel$lambdacoef,nclass),nrow=nclass,byrow=TRUE)
          initmodel <- fitAdaptRandom(patterns,freq=freq,
                                                  nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
                                                  initclassp=initmodel$classp,
                                                  initlambdacoef=initlambdacoef,
                                                  gh=norm.gauss.hermite(quadpoints),
                                                  constload=constload,blocksize=level2size,
                                                  probit=probit,byclass=byclass,qniterations=qniterations,
                                                  penalty=penalty,verbose=verbose)
        }
        # now fit the level 2
        initmodel <- fitAdaptRandom2(patterns,freq=freq,
                                                nclass=nclass,calcSE=calcSE,initoutcomep=initmodel$outcomep,
                                                initclassp=initmodel$classp,
                                                initlambdacoef=initmodel$lambdacoef,
                                                 initltaucoef=NULL,
                                                gh=norm.gauss.hermite(quadpoints),
                                                constload=constload,level2size=level2size,
                                                probit=probit,byclass=byclass,qniterations=qniterations,
                                                penalty=penalty,verbose=verbose)
      }
    }
# check rank of Hessian
    if (calcSE) {
      if (rankMatrix(initmodel$fit$hessian) < dim(initmodel$fit$hessian)[1])
        warning("Rank of Hessian less than number of estimated parameters - model is possibly underidentified or a parameter estimate is on the boundary")
    }
    fit <- initmodel
    fit$call <- cl
    fit$nclass <- nclass
    fit$random <- random
    fit$constload <- constload
    fit$level2 <- level2
    fit$level2size <- level2size
    fit$byclass <- byclass
    fit$probit <- probit
    fit$quadpoints <- quadpoints
    fit$blocksize <- blocksize
     fit$patterns <- patterns
    fit$notrials <- notrials
    fit$freq <- freq
    fit$qniterations <- qniterations
    fit$penalty <- penalty
    class(fit) <- "randomLCA"
    return(fit)
  }

