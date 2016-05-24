fitme <- function(formula,data, ## matches minimal call of HLfit
                             init=list(),
                             init.HLfit=list(),
                             ranFix=list(), 
                             lower=list(),upper=list(),
                             resid.model=~1,
                             control.dist=list(),
                             control=list(), ## optim.scale, nloptr, maxcorners
                             processed=NULL, 
                             family=gaussian(),
                             ...
) {
  oricall <- mc <- match.call(expand.dots=TRUE) ## mc including dotlist
  
  if (FALSE) { ## 
    mc$init.corrHLfit <- mc$init
    mc$init <- NULL
    mc$control.corrHLfit <- mc$control
    mc$control <- NULL
    ## mc may include non-default 'objective' and 'trace' for corrHLfit 
    mc[[1L]] <- quote(spaMM::corrHLfit)
    hlcor <- eval(mc,sys.frame(sys.parent()))
    return(hlcor)
  }
  
  ## Preventing confusions
  if (!is.null(mc$ranPars)) {
    stop("incorrect 'ranPars' argument in fitme() call. Use ranFix (ranPars is for HLCor() only)")
  }
  if ( ! (is.list(lower) && is.list(upper))) {
    wrongclass <- setdiff(unique(c(class(lower),class(upper))),"list")
    stop(paste("'lower' and 'upper' must be of class list, not",paste(wrongclass,collapse=" or ")))
  } ## as.list() would flatten rho vectors
  #
  HLnames <- (c(names(formals(HLCor)),names(formals(HLfit)),
                names(formals(designL.from.Corr)),names(formals(make_scaled_dist))))  
  dotnames <- setdiff(names(mc)[-1],names(formals(fitme)))
  argcheck <- setdiff(dotnames,HLnames)
  if (length(argcheck)>0) warning(paste("suspect argument(s) ",paste(argcheck, collapse=",")," in fitme call."))
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
  
  mc[[1L]] <- quote(spaMM::fitme_body) 
  hlcor <- eval(mc,parent.frame()) 
  attr(hlcor,"corrHLfitcall") <- oricall ## this says the hlcor was returned by corrHLfit (which should be fitme?)
  return(hlcor)
}



fitme_body <- function(processed,
                       init=list(),
                       init.HLfit=list(),
                       ranFix=list(), 
                       lower=list(),upper=list(),
                       control.dist=list(),
                       control=list(), ## optim.scale, Optimizer, optimizer.args, maxIter, maxcorners, precision
                       verbose, ## provided by corrHLfit
                       ... ## cf dotnames processing below
) {
  
  dotlist <- list(...) ## forces evaluations, which makes programming easier...
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
  # 'processed' may be modified below, then will be copied in HLCor.args (and then removed from this envir for safety)
  HLCor.args$verbose <- verbose ## fitme_body fn argument
  
  RHOMAX <- NULL
  NUMAX <- 50
  
  init.optim <- init ## usages of init$rho, etc. must be tracked through init.optim too  
  optim.scale <- control$optim.scale 
  if (is.null(optim.scale)) optim.scale="transformed" ## currently no public alternative
  #
  maxcorners <- control$maxcorners
  ## 2015/12/24 : setting maxcorners to 0L till makes a difference for
  # test-Nugget -> different p_v / p_bv compromise (!)
  # nullfit of blackcap (test-fixedLRT) is better with the corners.
  if (is.null(maxcorners)) maxcorners <- 2^11
  #
  spatial.terms <- findSpatial(getProcessed(processed,"predictor",from=1L))
  spatial.model <- spatial.terms[[1L]] 
  corr.model <- paste("",as.character(spatial.model[[1L]]),sep="") ## so that LHS= "" if RHS= spatial.model[[1L]]
  HLCor.args$corr.model <- corr.model 

  ## distMatrix or uniqueGeo potentially added to HLCor.args:
  ############### (almost) always check geo info ###################
  # modify HLCor.args and <>bounds
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
    if ( inherits(data,"list")) {
      dataForCheck <- data[[1]]
    } else dataForCheck <- data
    if ( ! is.null(spatial.model)) coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(dataForCheck))
  } 
  # 
  Fixrho <- getPar(ranFix,"rho")
  if ( (! is.null(Fixrho)) && (! is.null(init$rho)) ) {
    stop("(!) 'rho' given as element of both 'ranFix' and 'init'. Check call.")    
  } else rho.size <- max(length(Fixrho),length(init$rho))
  if ( (! is.null(getPar(ranFix,"nu"))) && (! is.null(init$nu)) ) {
    stop("(!) 'nu' given as element of both 'ranFix' and 'init'. Check call.")    
  }
  if ( (! is.null(getPar(ranFix,"ARphi"))) && (! is.null(init$ARphi)) ) {
    stop("(!) 'ARphi' given as element of both 'ranFix' and 'init'. Check call.")    
  }
  if ( (! is.null(getPar(ranFix,"Nugget"))) && (! is.null(init$Nugget)) ) {
    stop("(!) 'Nugget' given as element of both 'ranFix' and 'init'. Check call.")    
  }
  #
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
  #
  init.optim$lambda <- 1 ## FR->FR DEVEL !!!!
  cat("fitme: init.optim$lambda <- 1 !")
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
  needHLCor_specific_args <- (length(setdiff(names(ranPars),c("trPhi","trLambda")))>0L) 
  if (needHLCor_specific_args) {
    HLcallfn.obj <- "HLCor.obj"
    HLcallfn <- "HLCor"
    HLCor.args$control.dist <- control.dist ## modified locally
  } else {
    HLcallfn.obj <- "HLfit.obj"
    HLcallfn <- "HLfit"
  }
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
  } else optimMethod <- "nloptr"
  HLCor.args$processed <- processed
  processed <- "'processed' erased after copy in 'HLCor.args' to make sure it is not modified later"
  ## 
  anyHLCor_obj_args <- HLCor.args
  ## HLCor.obj uses a vector + skeleton
  if (needHLCor_specific_args) {
    anyHLCor_obj_args$skeleton <- structure(init.optim,RHOMAX=RHOMAX,NUMAX=NUMAX) ## logscale, only used by HLCor.obj
  } else anyHLCor_obj_args$skeleton <- init.optim
  attr(anyHLCor_obj_args$skeleton,"type") <- list() ## declares a list of typeS of elemnts of skeleton
  attr(anyHLCor_obj_args$skeleton,"type")[names(init.optim)] <- "fix" # fixed with the HLCor call 
  if (ncol(processed$X.Re)>0) {
    anyHLCor_obj_args$`HLCor.obj.value` <- "p_bv" 
  } else anyHLCor_obj_args$`HLCor.obj.value` <- "p_v" 
  initvec <- unlist(init.optim)
  ####    tmpName <- generateName("HLtmp") ## tmpName is a string such as "HLtmp0"
  #    anyHLCor_obj_args$init.HLfit <- tmpName 
  ####    assign(tmpName,list(),pos=".GlobalEnv") ## sets HLtmp0 (or a similarly named variable) at the global level
  if (optimMethod=="iterateSEMSmooth") {
    stop("reimplement later")
  } else { ## this is also called if length(lower)=0 by  (SEM or not) and optPars is then null 
    objfn_nloptr <- function(x,anyHLCor_obj_args) { ## all functions should have the same args.
      arglist <- c(list(ranefParsVec=x),anyHLCor_obj_args)
      return( - do.call(HLcallfn.obj,arglist))
    }
    nloptr_controls <- list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1.0e-4,maxeval=-1,print_level=0)
    nloptr_controls[names(control$nloptr)] <- control$nloptr
    optr <- nloptr(x0=unlist(init.optim),eval_f=objfn_nloptr,lb=unlist(lower),ub=unlist(upper),
                   opts=nloptr_controls,anyHLCor_obj_args=anyHLCor_obj_args)
    optPars <- structure(optr$solution,method="nloptr",optr=optr)
  }
  ranPars[names(init.optim)] <- optPars ## avoids overwriting fixed ranPars 
  attr(ranPars,"type")[names(init.optim)] <- "outer" ##  
  if (needHLCor_specific_args) {
    attr(ranPars,"RHOMAX") <- RHOMAX
    attr(ranPars,"NUMAX") <- NUMAX
  }
  HLCor.args$ranPars <- ranPars ## variable locally
  verbose["warn"] <- TRUE ## important!
  HLCor.args$verbose <- verbose ## modified locally
  hlcor <- do.call(HLcallfn,HLCor.args) ## recomputation post optimization (or only computation, if length(lower)=0)
  if (needHLCor_specific_args) {
    if ( is.null(HLCor.args$adjMatrix) && is.null(attr(hlcor,"info.uniqueGeo")) ) { ## typically if DistMatrix was passed to HLCor...
      attr(hlcor,"info.uniqueGeo") <- uniqueGeo ## uniqueGeo should have been computed in all relevant cases where test is true (tricky)
    }
  }
  attr(hlcor,"objective") <- anyHLCor_obj_args$`HLCor.obj.value` 
  ## 
  optimInfo <- list(optim.pars=optr$par, 
                                  init.optim=init.optim,
                                  lower=lower,upper=upper,
                                  user.lower=user.lower,user.upper=user.upper) 
  if (needHLCor_specific_args) {optimInfo$RHOMAX <- RHOMAX}
  attr(hlcor,"optimInfo") <- optimInfo ## RHOMAX is also in the *call* of the hlcor object...
  return(hlcor)
}