## wrapper for optimization of HLCor.obj OR (iterateSEMSmooth -> HLCor directly)
corrHLfit <- function(formula,data, ## matches minimal call of HLfit
                      init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      trace=list(file=NULL,append=TRUE),
                      objective="p_bv", ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      resid.model=~1,resid.formula,
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, maxcorners, precision
                      processed=NULL, ## added 2014/02 for programming purposes
                      family=gaussian(),
                      ... ## pb est risque de passer des args mvs genre HL.method et non HLmethod...
) {
  oricall <- mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  ## Preventing confusions
  if (!is.null(mc$ranPars)) {
    stop("incorrect 'ranPars' argument in corrHLfit call. Use ranFix (ranPars is for HLCor only)")
  }
  if ( ! (is.list(lower) && is.list(upper))) {
    wrongclass <- setdiff(unique(c(class(lower),class(upper))),"list")
    stop(paste("'lower' and 'upper' must be of class list, not",paste(wrongclass,collapse=" or ")))
  } ## as.list() would flatten rho vectors
  #
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(corrHLfit)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit call."))
  # 
  if (is.null(processed)) {
    family <- checkRespFam(family)
    FHF <- formals(HLfit) ## makes sure about default values 
    names_FHF <- names(FHF)
    if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
    names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
    FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
    preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
    preprocess.formal.args$family <- family ## already checked 
    preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
    preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
    #
    famfam <- family$family
    if ( identical(family$family,"multi")) {
      ## then data are reformatted as a list. Both HLCor and HLfit can analyse such lists for given corrPars and return the joint likelihood
      ## By contrast HLCor should not fit different corrPars to each data, so it does not lapply("corrHLfit",...)
      ## Rather, it calls preprocess which will construct a list of processed objects, to be used conjointly with the data list.
      ## But then we should not attempt to modify an element of 'pocessed' as if it was a single processed object
      ## We must use setProcessed / getProcessed to access list elements.
      if ( ! inherits(data,"list")) {
        familyargs <- family
        familyargs$family <- NULL
        familyargs$binfamily <- NULL
        ## we need the data list in the corrHLfit envir for the call to makeCheckGeoMatrices
        preprocess.formal.args$data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
      }     
    }
    mc$processed <- do.call(preprocess,preprocess.formal.args,envir=environment(formula))
    mc$verbose <- reformat_verbose(mc$verbose,For="corrHLfit")
    ## removing all elements that are matched in processed:
    mc$data <- NULL
    mc$family <- NULL
    mc$formula <- NULL ## processed
    mc$HLmethod <- NULL ## processed$HL  
  }  
  
  mc[[1L]] <- quote(spaMM::corrHLfit_body) 
  hlcor <- eval(mc,parent.frame()) 
  attr(hlcor,"corrHLfitcall") <- oricall ## this says the hlcor was returned by corrHLfit
  return(hlcor)
}


corrHLfit_body <- function(processed,
                           init.corrHLfit=list(),
                      init.HLfit=list(),
                      ranFix=list(), 
                      lower=list(),upper=list(),
                      trace=list(file=NULL,append=TRUE),
                      objective="p_bv", ## return value of HLCor.obj for optim calls... FR->FR meaningless for full SEM
                      control.dist=list(),
                      control.corrHLfit=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, maxcorners, precision
                      verbose, ## provided by corrHLfit
                      ... ## cf dotnames processing below
) { 

  

  #########################################################
  #
  dotlist <- list(...) ## forces evaluations, which makes programming easier...

  if (FALSE) { ## devel code
    HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                  names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
    argcheck <- setdiff(names(dotlist),HLnames)
    if (length(argcheck)>0) stop(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in corrHLfit_body call."))
  }
  
  data <- getProcessed(processed,"data") ## gets a data list if processed is a list
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  ## cf parallel code in HLCor.obj
  ## fill HLCor.args
  good_dotnames <- intersect(names(dotlist),HLnames) ## those specifically for the called fns as def'd by HLnames
  if (length(good_dotnames)>0L) {
    HLCor.args <- dotlist[good_dotnames]
  } else HLCor.args <- list() 
  typelist <- list() ## on veut une list pour pouvoir supp des elements par <- NULL
  typelist[names(ranFix)] <- "fix"
  attr(ranFix,"type") <- typelist
  ## replace some HLCor.args members  
  if ( ! is.null(attr(processed,"multiple"))) {
    pnames <- names(processed[[1]])
  } else pnames <- names(processed)
  for (st in pnames) HLCor.args[st] <- NULL 
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety) => NO NEED to copy any processed element in HLCor.args now.
  HLCor.args$verbose <- verbose ## corrHLfit_body fn argument
  
  RHOMAX <- NULL
  NUMAX <- 50
  
  init.optim <- init.corrHLfit ## usages of init.corrHLfit$rho, etc. must be tracked through init.optim too  
  optim.scale <- control.corrHLfit$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  Optimizer <- control.corrHLfit$Optimizer ## either "nlminb" or one of the methods of optim()
  if (is.null(Optimizer)) Optimizer="L-BFGS-B" ## default for locoptim but still requested for "Optimizer=Optimizer"
  optimizers.args <- control.corrHLfit[c("nlminb","optim","optimize")] 
  maxcorners <- control.corrHLfit$maxcorners
  ## 2015/12/24 : setting maxcorners to 0L till makes a difference for
  # test-Nugget -> different p_v / p_bv compromise (!)
  # nullfit of blackcap (test-fixedLRT) is better with the corners.
  if (is.null(maxcorners)) maxcorners <- 2^11 
  maxIter <- control.corrHLfit$maxIter 
  if (is.null(maxIter)) maxIter<- 10000
  if ( ! (objective %in% c("p_v","p_bv"))) {
    mess <- pastefrom("invalid value of the 'objective' argument.",prefix="(!) From ")
    stop(mess)
  }
  if (length(trace)>0 && ! all(names(trace) %in% c("file","append"))) {
    mess <- pastefrom("'trace' elements should be named and = 'file','append'.",prefix="(!) From ")
    message(mess)
  } 
  if (is.character(trace$file) && ( ! trace$append) ) {
    try(unlink(trace$file))
  }  ## the file is written in by HLCor()                 
  HLCor.args$trace <- trace$file
  #
  spatial.terms <- findSpatial(getProcessed(processed,"predictor",from=1L))
  spatial.model <- spatial.terms[[1]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1]]) 
  } else {
    # corr.model <- "Matern" ## up to v 1.0; for the defective syntax of the scripts for the Ecography paper
    stop("spatial correlation model not specified in 'formula': was valid in version 1.0 but not later.")
  } 
  HLCor.args$corr.model <- corr.model ## determined locally
  #
  ## distMatrix or uniqueGeo potentially added to HLCor.args:
  ############### (almost) always check geo info ###################
  if (corr.model  %in% c("SAR_WWt","adjacency","ar1")) {
    ## adjMatrix should become weightMatrix
    if ( is.null(HLCor.args$adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
    if (isSymmetric(HLCor.args$adjMatrix)) {
      decomp <- selfAdjointWrapper(HLCor.args$adjMatrix)
      attr(HLCor.args$adjMatrix,"symSVD") <- decomp
    } else {
      if (corr.model  %in% c("SAR_WWt")) {
        decomp <- eigen(HLCor.args$adjMatrix,symmetric=FALSE) ## FR->FR not RcppEigen-optimized
        attr(HLCor.args$adjMatrix,"UDU.") <- list(u=decomp$vectors,d=decomp$values,u.=solve(decomp$vectors))
      } else stop("'adjMatrix' is not symmetric") ## => invalid cov mat for MVN
    }
    rhorange <- sort(1/range(decomp$d)) ## keeping in mind that the bounds can be <>0
    if(verbose["SEM"]) cat(paste("Feasible rho range: ",paste(signif(rhorange,6),collapse=" -- "),"\n"))
  } else {
    if ( is.null(spatial.model)) {
      stop("An obsolete syntax for the adjacency model appears to be used.")
      ## coordinates <- c("x","y") ## backward compatibility... old syntax with (1|pos) and default values of the coordinates argument
    } else {
      if ( inherits(data,"list")) {
        dataForCheck <- data[[1]]
      } else dataForCheck <- data
      coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(dataForCheck))
    }
  } ## 
  #
  Fixrho <- getPar(ranFix,"rho")
  if ( (! is.null(Fixrho)) && (! is.null(init.corrHLfit$rho)) ) {
    stop("(!) 'rho' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  } else rho.size <- max(length(Fixrho),length(init.corrHLfit$rho))
  if ( (! is.null(getPar(ranFix,"nu"))) && (! is.null(init.corrHLfit$nu)) ) {
    stop("(!) 'nu' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if ( (! is.null(getPar(ranFix,"ARphi"))) && (! is.null(init.corrHLfit$ARphi)) ) {
    stop("(!) 'ARphi' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if ( (! is.null(getPar(ranFix,"Nugget"))) && (! is.null(init.corrHLfit$Nugget)) ) {
    stop("(!) 'Nugget' given as element of both 'ranFix' and 'init.corrHLfit'. Check call.")    
  }
  if (corr.model %in% c("Matern","AR1")) {
    rho_mapping <- control.dist$rho.mapping ## may be NULL
    if (is.null(rho_mapping) ) { ## && length(coordinates)>1L ?
      if (length(coordinates)==rho.size) { ## rho.size comes from explicit rho from user
        rho_mapping <- seq_len(rho.size)           
        names(rho_mapping) <- coordinates
        control.dist$rho.mapping <- rho_mapping
      } else if (length(rho.size)>1L) stop("'rho.mapping' missing with no obvious default from the other arguments.")
    } ## then (for given corr.model's) there is rho_mapping
    locarglist<- list(data=data,distMatrix=dotlist$distMatrix,
                      uniqueGeo=HLCor.args$uniqueGeo,coordinates=coordinates)
    if(!is.null(dist.method <- control.dist$dist.method)) locarglist$dist.method <- dist.method
    geoMats <- do.call(makeCheckGeoMatrices,locarglist)
    nbUnique <- geoMats$nbUnique  
    if (inherits(nbUnique,"list")) nbUnique <- mean(unlist(nbUnique))
    distMatrix <- geoMats$distMatrix   
    uniqueGeo <- geoMats$uniqueGeo   
    if (rho.size<2) { ## can be 0 if no explicit rho in the input  
      HLCor.args$distMatrix <- distMatrix   ## determined locally
    } else {
      HLCor.args$uniqueGeo <- uniqueGeo   ## determined locally
    }
  } else nbUnique <- NULL
  
  #
  arglist <- list(init.optim=init.optim,init.HLfit=init.HLfit,ranFix=ranFix,
                  corr.model=corr.model,optim.scale=optim.scale,
                  control.dist=control.dist)
  if (corr.model %in% c("Matern")) {
    maxrange <- calc_maxrange(rho.size,distMatrix,uniqueGeo,rho_mapping,dist.method) 
    RHOMAX <- 1.1*30*nbUnique/maxrange ## matches init$rho in calc_inits()
    arglist <-c(arglist,list(maxrange=maxrange,RHOMAX=RHOMAX,NUMAX=NUMAX))
  }
  if (corr.model %in% c("SAR_WWt","adjacency","ar1") &&  is.null(getPar(ranFix,"rho"))
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
  ) arglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  inits <- do.call("calc_inits",arglist)
  init <- inits$`init` ## keeps all init values, all in untransformed scale
  init.optim <- inits$`init.optim` ## subset of all optim estimands, as name implies, and in transformed scale
  init.HLfit <- inits$`init.HLfit` ## subset as name implies 
  #
  
  ## maximization VIA HLCor.obj
  ## by maxim over (corrpars,phi,lambda, (beta)_[corrpars,phi,lambda])
  ##     if trPhi,trLambda are in the init.optims
  ## or by maxim over (corrpars,(beta,phi,lambda)_corrpars)
  ##     otherwise.
  ################ construct intervals for this maximization
  ## construct default upper and lower values ; on transformed scale by default
  user.lower <- lower; user.upper <- upper ## keep user input 
  if ("lambda" %in% c(names(user.lower),names(user.lower)) 
      && is.null(init$lambda)) {
    stop("'lambda' in 'lower' or 'upper' has no effect if absent from 'init.corrHLfit'.")
  }
  ################
  LUarglist <- list(canon.init=init,
                    init.optim=init.optim, ## initially with right transformed variables but wrong values
                    user.lower=user.lower,user.upper=user.upper,
                    corr.model=corr.model,nbUnique=nbUnique,
                    ranFix=ranFix,control.dist=control.dist,
                    optim.scale=optim.scale, RHOMAX=RHOMAX,NUMAX=NUMAX)
  if (corr.model %in% c("SAR_WWt","adjacency","ar1") &&  is.null(getPar(ranFix,"rho"))
      && (! is.numeric(init.HLfit$rho)) ## init.HLfit$rho$rho NULL or NA) 
      ) { 
    LUarglist$rhorange <- rhorange ## will serve to initialize either HLfit or optim 
  }
  LowUp <- do.call("makeLowerUpper",LUarglist)
  ## LowUp: a list with elements lower and upper that inherits names from init.optim, must be optim.scale as init.optim is by construction
  lower <- LowUp$lower ## list ! which elements may have length >1 !
  upper <- LowUp$upper ## list !
  #
  ranPars <- ranFix ## ranPars argument of HLCor contains both fixed and estimated parameters:
  varNames <- names(init.HLfit) ## hence those that will be variable within HLfit
  ranPars[varNames] <- init.HLfit[varNames] ## FR->FR duplicat (?) qui montre qu'un attribute serait mieux
  attr(ranPars,"type")[varNames] <- "var"  
  HLCor.args$ranPars <- ranPars  ## variable locally
  #
  HLCor.args$control.dist <- control.dist ## modified locally
  processedHL1 <- getProcessed(processed,"HL[1]",from=1L) ## there's also HLmethod in processed<[[]]>$callargs
  if (!is.null(processedHL1) && processedHL1=="SEM" && length(lower)>0) {
    optimMethod <- "iterateSEMSmooth"
    processed <- setProcessed(processed,"SEMargs$SEMseed",value="NULL") ## removing default SEMseed
    ## : bc SEMseed OK to control individual SEMs but not  series of SEM 
    if (is.null(getProcessed(processed,"SEMargs$control_pmvnorm$maxpts",from=1L))) {
      if (length(LowUp$lower>0L)) {
        processed <- setProcessed(processed,"SEMargs$control_pmvnorm$maxpts",value="quote(250L*nrow(pmvnorm.Sig))") ## removing default SEMseed
      } ## else default visible in SEMbetalambda
    }
  } else optimMethod <- "locoptim"
  HLCor.args$processed <- processed
  processed <- "'processed' erased after copy in 'HLCor.args' to make sure it is not modified later"
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  anyHLCor_obj_args$skeleton <- structure(init.optim,RHOMAX=RHOMAX,NUMAX=NUMAX) ## logscale, only used by HLCor.obj
  attr(anyHLCor_obj_args$skeleton,"type") <- list() ## declares a list of typeS of elemnts of skeleton
  attr(anyHLCor_obj_args$skeleton,"type")[names(init.optim)] <- "fix" # fixed with the HLCor call 
  anyHLCor_obj_args$`HLCor.obj.value` <- objective ## p_v when fixedLRT-> corrHLfit
  anyHLCor_obj_args$traceFileName <- trace$file
  initvec <- unlist(init.optim)
  ####    tmpName <- generateName("HLtmp") ## tmpName is a string such as "HLtmp0"
  #    anyHLCor_obj_args$init.HLfit <- tmpName 
  ####    assign(tmpName,list(),pos=".GlobalEnv") ## sets HLtmp0 (or a similarly named variable) at the global level
  if (optimMethod=="iterateSEMSmooth") {   
    MAX <- list(trRho=RHOMAX, trNu=NUMAX) ##  MAX is used in the SEMdiagnosticplot...;
    ## its names should match the colnames of the data in Krigobj = the  parameters of the likelihood surface. Current code maybe not general.
    loclist <- list(anyHLCor_obj_args=anyHLCor_obj_args,  ## contains $processed
                       LowUp=LowUp,init.corrHLfit=init.corrHLfit, 
                       #preprocess.formal.args=preprocess.formal.args, 
                       control.corrHLfit=control.corrHLfit,
                       MAX=MAX)
    optr <- do.call("iterateSEMSmooth",loclist)
    optPars <- as.list(optr$par)
    if (!is.null(optPars)) attr(optPars,"method") <-"optimthroughSmooth"
  } else if (FALSE) { ## renewed coding of the iterative algo (only p_v); not documented => not checked for a long time
    optPars <- alternating(init.optim=init.optim,LowUp=LowUp,
                           anyHLCor_obj_args=anyHLCor_obj_args,maxIter=maxIter,
                           ranPars=ranPars,HLCor.args=HLCor.args,trace=trace,Optimizer=Optimizer,
                           optimizers.args=optimizers.args,maxcorners=maxcorners)
    if (!is.null(optPars)) attr(optPars,"method") <-"alternating"
  } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null 
    loclist<-list(init.optim=init.optim,LowUp=LowUp,anyObjfnCall.args=anyHLCor_obj_args,trace=trace,Optimizer=Optimizer,
                  optimizers.args=optimizers.args,maxcorners=maxcorners,maximize=TRUE) 
    optPars <- do.call("locoptim",loclist)
    if (!is.null(optPars)) attr(optPars,"method") <- "locoptim"
  }
  ranPars[names(optPars)] <- optPars ## avoids overwriting fixed ranPars 
  attr(ranPars,"type")[names(optPars)] <- "outer" ##  
  attr(ranPars,"RHOMAX") <- RHOMAX
  attr(ranPars,"NUMAX") <- NUMAX
  HLCor.args$ranPars <- ranPars ## variable locally
  verbose["warn"] <- TRUE ## important!
  HLCor.args$verbose <- verbose ## modified locally
  hlcor <- do.call("HLCor",HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  if ( is.null(HLCor.args$adjMatrix) && is.null(attr(hlcor,"info.uniqueGeo")) ) { ## typically if DistMatrix was passed to HLCor...
    attr(hlcor,"info.uniqueGeo") <- uniqueGeo ## uniqueGeo should have been computed in all relevant cases where test is true (tricky)
  }
  attr(hlcor,"objective") <- anyHLCor_obj_args$`HLCor.obj.value` 
  ## 
  attr(hlcor,"optimInfo") <- list(optim.pars=optPars, 
                                  init.optim=init.optim,
                                  lower=lower,upper=upper,
                                  user.lower=user.lower,user.upper=user.upper,
                                  RHOMAX=RHOMAX) ## RHOMAX is also in the *call* of the hlcor object...
  if ( ! is.null(locoptr <- attr(optPars,"optr")) && locoptr$convergence>0L) {
    hlcor$warnings$optimMessage <- paste("optim() message: ",locoptr$message," (convergence=",locoptr$convergence,")",sep="")
  }
  if ( ( ! is.null(optPars)) && attr(optPars,"method")== "optimthroughSmooth") {
    # provide logL estimate from the smoothing, to be used rather than the hlcor logL :
    logLapp <- optr$value
    attr(logLapp,"method") <- "  logL (smoothed)" 
    if (FALSE) { ## essai correction biais
      message("Trying to correct bias of SEM procedure...")
      p_lambda <- NROW(hlcor$lambda.object$coefficients_lambda)
      np <- NCOL(hlcor$`X.pv`) + p_lambda
      lmframe <- replicate(4*((np+1)*(np+2)/2), {## ((np+1)*(np+2)/2) is number of params to be fitted
        rephlcor <- do.call("HLCor",HLCor.args)
        if (p_lambda==0L) {
          c(fixef(rephlcor), logLik(rephlcor)[1])
        } else if (p_lambda==1L) {
          c(fixef(rephlcor), lambda=rephlcor$lambda, logLik(rephlcor)[1])
        } else stop("code missing for p_lambda>0 for bias correction of SEM.")
        ## FR->FR + problem if one of the predict vaiables in named logLapp
      } )
      lmframe <-as.data.frame(t(lmframe))
      colnames(lmframe)[colnames(lmframe)=="(Intercept)"] <- "fixef.intercept" ## because "(Intercept)" does not work below
      quadfit <- eval(parse(text=paste("lm(logLapp ~ polym(",
                                       paste(colnames(lmframe)[1:np],collapse=","),
                                       ", degree=2,raw=TRUE),data=lmframe)"))) ## predict does not work if raw=FALSE
      lmlower <- apply(lmframe[,1:np,drop=FALSE],2,min)
      lmupper <- apply(lmframe[,1:np,drop=FALSE],2,max)
      parnames <- names(lmlower)
      ## R issue: [.data.frame loses col name for single col.
      ## => lmframe[which.max(lmframe$logLapp),1:np] is 1-row data.frame is np>1; is *unnamed* numeric if np=1 unless drop=FALSE; 
      ##    where we would wish it to be named numeric in both cases
      initpar <- unlist(lmframe[which.max(lmframe$logLapp),1:np,drop=FALSE]) 
      optpolym <- optim(par=initpar, 
                        fn=function(v) {
                          v <- data.frame(matrix(v,nrow=1)) ## matrix() loses names... ## data.frame uses names as col or row names dependeing on length(v)...
                          colnames(v) <- parnames ## hence give names as a last step; again drop=F will be needed for 1-col 
                          predict(quadfit, newdata=v[c(1,1),,drop=FALSE])[1] ## [c(1,1),] is an awful patch for predict.poly
                        }, 
                        lower=lmlower,upper=lmupper,method="L-BFGS-B",
                        control=list(fnscale=-1,parscale=lmupper-lmlower)) 
      ## computation of dfs as in compare.model.structures
      attr(logLapp,"optpolym") <- optpolym$value
      attr(logLapp,"lmframe") <- lmframe
      attr(logLapp,"replicates") <- replicate(20,logLik(do.call("HLCor",HLCor.args)))
      attr(logLapp,"stddevres") <- ((predict(optr$Krigobj)-optr$Krigobj$data$logLobj)/(1-optr$Krigobj$lev_phi))
      attr(logLapp,"df_replicates") <- np
    }
    hlcor$APHLs$logLapp <- logLapp
  }
  if (is.character(trace$file)) {
    ## crude display of variable names in the trace file
    traceNames <- paste("# ",paste(names(hlcor$APHLs),collapse=" "))
    traceNames <- paste(traceNames,"lambda",sep=" ")
    if ( ! is.null(hlcor$phi)) traceNames <- paste(traceNames,"phi",sep=" ")
    traceNames <- paste(traceNames,paste(names(anyHLCor_obj_args$skeleton),collapse=" "),sep=" ")
    traceNames <- paste(traceNames," and optim parameters in canonical scale  ",sep=" ")
    write(traceNames,file=trace$file,append=T)   
  }
  ####  rm(list=c(tmpName),pos=".GlobalEnv") ## removes HLtmp0 at the global level
  return(hlcor) ## it's the call which says it was returned by corrHLfit
}

    
