## the following must match the'unique' method is ULI as explained there
calcUniqueGeo <- function(data) {
  redondGeo <- apply(data,1,paste,collapse=" ") ## creates character string
  dfforunique <- cbind(data,redondGeo) ## associates rownames of data to redondGeo
  uniqueGeo <- unique(dfforunique[,ncol(dfforunique),drop=FALSE]) ## keeps rownames of first instances
  uniqueGeo <- data[rownames(uniqueGeo),,drop=FALSE] ## uses rownames, 'unique' numeric values based on character representations 
  return(uniqueGeo)
}


`extract.check.coords` <- function(spatial.model,datanames) {
  if ( ! is.null(spatial.model)) {
    bars <- spatial.model[[2]] 
    coordinates <- DEPARSE(bars[[3]]) ## "x + y"
    coordinates <-  strsplit(coordinates," ")[[1]]
    coordinates <- setdiff(coordinates,c("+","%in%",":","/","")) ## "" for hidden linebreaks (?)
  } else {
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
    ## very old code handling old syntax with (1|pos) and default values of the coordinates argument
    coordinates <- c("x","y") ## back compat
  }
  coordcols <- which(datanames %in% coordinates)
  if ( length(coordcols) != length(coordinates) ) {
    stop("variables 'coordinates' not all found in the 'data'")
  }
  return(coordinates) ## should be ordered as bars[[3]] (important for predict)
}

## better for development to avoid name conflicts with OKsmooth :toCanonical and :canonize
canonizeRanPars <- function(ranPars, ## should have a RHOMAX attribute when trRho in input
                            corr.model,checkComplete=TRUE) {
  trueCorrpars <- list()
  if (corr.model %in% c("Matern")) {
    if (!is.null(ranPars$trNu)) { ## either we have nu,rho or trNu,trRho 
      ranPars$nu <- nuInv(ranPars$trNu,ranPars$trRho,NUMAX=attr(ranPars,"NUMAX")) ## before trRho is removed...
      ranPars$trNu <- NULL
      attr(ranPars,"type")$nu <- attr(ranPars,"type")$trNu
      attr(ranPars,"type")$trNu <- NULL
    } 
    nu <- ranPars$nu
    if (is.null(nu) && checkComplete) {
      mess <- pastefrom("nu missing from ranPars (or correlation model mis-identified).",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$nu <- nu 
  } 
  if (corr.model=="AR1") {
    ARphi <- ranPars$ARphi
    if (is.null(ARphi) && checkComplete) {
      mess <- pastefrom("ARphi missing from ranPars.",prefix="(!) From ")
      stop(mess)
    }
    trueCorrpars$ARphi <- ARphi    
  } else if (corr.model != "corrMatrix") { ## all models with a 'rho' parameter
    if (!is.null(ranPars$trRho)) {
      ranPars$rho <- rhoInv(ranPars$trRho,RHOMAX=attr(ranPars,"RHOMAX"))  
      ranPars$trRho <- NULL
      attr(ranPars,"type")$rho <- attr(ranPars,"type")$trRho
      attr(ranPars,"type")$trRho <- NULL
    } ## else there may simply be rho rather than trRho (including for adjacency model through optim procedure !)
    trueCorrpars$rho <- rho <- ranPars$rho
    if (is.null(rho)) {
      if(corr.model=="adjacency") { ## then allow a direct call through HLCor 
        ranPars$rho <- 0
        attr(ranPars,"type")$rho <- "var"
      } else if (checkComplete) {
        mess <- pastefrom("rho missing from ranPars.",prefix="(!) From ")
        stop(mess)
      }
    } 
  }
  Nugget <- ranPars$Nugget
  if (! is.null(Nugget)) trueCorrpars$Nugget <- Nugget 
  if (!is.null(ranPars$trPhi)) {
    ranPars$phi <- dispInv(ranPars$trPhi)
    ranPars$trPhi <- NULL
    attr(ranPars,"type")$phi <- attr(ranPars,"type")$trPhi
    attr(ranPars,"type")$trPhi <- NULL
  } else if (!is.null(ranPars$logphi)) { ## debug code
    ## HL.info$ranFix$phi <- exp(ranPars$logphi)
    stop("logphi in HLCor...")
  } #####################  else HL.info$ranFix$phi <- ranPars$phi ## y st deja !?
  if (!is.null(ranPars$trLambda)) {## 
    ranPars$lambda <- dispInv(ranPars$trLambda)
    ranPars$trLambda <- NULL
    attr(ranPars,"type")$lambda <- attr(ranPars,"type")$trLambda
    attr(ranPars,"type")$trLambda <- NULL
  } else if (!is.null(ranPars$loglambda)) { ## debug code
    stop("loglambda in HLCor...")
  } ##################### else HL.info$ranFix$lambda <- ranPars$lambda
  return(list(trueCorrpars=trueCorrpars,ranPars=ranPars))
}


HLCor <- function(formula,
                  data,family=gaussian(),
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),control.dist=list(),
                  ...) { 
  oricall <- mc <- match.call(expand.dots = TRUE) 
  if ( ! is.null(mc$ranFix)) { ## avoiding user's confusion
    stop("!From HLCor: ranFix found in '...'. Make sure to use ranPars only")
  }
  if (!is.null(mc$LamFix)) {
    stop("argument LamFix of HLCor is obsolete")
  }
  if (!is.null(mc$PhiFix)) {
    stop("argument PhiFix of HLCor is obsolete")
  }  
  # frst steps as in HLFit: (no need to test missing(data) in several functions)
  if (is.null(processed <- mc$processed)) { ## no 'processed'
    ## FR->FR suggests we should add processed as argument of HLCor...
    family <- checkRespFam(family)
    if ( identical(family$family,"multi")) {
      if ( ! inherits(data,"list")) {
        if(family$binfamily$family=="binomial") {
          familyargs <- family
          familyargs$family <- NULL
          familyargs$binfamily <- NULL
          data <- do.call(binomialize,c(list(data=data),familyargs)) ## if data not already binomialized
        }
      }
    }
    if ( inherits(data,"list")) {
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(data)),function(it){
        locmc <- mc
        if (identical(family$family,"multi")) locmc$family <- family$binfamily
        locmc$data <- data[[it]]
        locmc$distMatrix <- mc$distMatrix[[it]]
        locmc$uniqueGeo <- mc$uniqueGeo[[it]]
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else {## there is a single data set, still without processed
      FHF <- formals(HLfit) ## makes sure about default values 
      names_FHF <- names(FHF)
      if ( ! is.null(mc$resid.formula)) mc$resid.model <- mc$resid.formula
      names_nondefault  <- intersect(names(mc),names_FHF) ## mc including dotlist
      FHF[names_nondefault] <- mc[names_nondefault] ##  full HLfit args
      preprocess.formal.args <- FHF[which(names_FHF %in% names(formals(preprocess)))] 
      preprocess.formal.args$family <- family ## already checked 
      preprocess.formal.args$rand.families <- FHF$rand.family ## because preprocess expects $rand.families 
      preprocess.formal.args$predictor <- FHF$formula ## because preprocess stll expects $predictor 
      preprocess.formal.args$ranFix <- ranPars ## because preprocess expects ranFix
      mc$processed <- do.call(preprocess,preprocess.formal.args,envir=environment(formula))
      # HLCor_body() called below
    }
  } else { ## 'processed' is available
    multiple <- attr(processed,"multiple")
    if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)),function(it){
        locmc <- mc
        locmc$processed <- processed[[it]] ## The data are in processed !
        locmc$distMatrix <- distMatrix[[it]] ## but the matrices are not HLfit args hence not in processed ! 
        locmc$uniqueGeo <- uniqueGeo[[it]]
        eval(locmc)
      }) ## a pure list of HLCor objects
      liks <- sapply(fitlist,function(v) {unlist(v$APHLs)})
      liks <- apply(liks,1,sum)
      attr(fitlist,"APHLs") <- as.list(liks)
      class(fitlist) <- c("HLfitlist",class(fitlist)) 
      return(fitlist) ## list of HLfit object + one attribute
    } else { ## there is one processed for a single data set 
      # HLCor_body() called below
    }
  }
  ################# single processed, single data analysis: 
  mc$verbose <- reformat_verbose(mc$verbose,For="HLCor")
  mc$data <- NULL
  mc$family <- NULL
  mc$formula <- NULL
  mc$prior.weights <- NULL
  mc$HLmethod <- NULL ## processed$HL  
  mc$rand.family <- NULL ## processed$rand.families  
  mc$resid.formula <- NULL ## mc$resid.model  
  mc$REMLformula <- NULL 
  mc[[1L]] <- quote(spaMM::HLCor_body)
  hlcor <- eval(mc,parent.frame())
  attr(hlcor,"HLCorcall") <- oricall ## potentially used by getCall(object) in update.HL 
  if (mc$verbose["HLCorSummary"]) { ## useful in final call from corrHLfit
    summary(hlcor) ## input corr pars have been printed at the beginning...   
  }
  return(hlcor)
}

HLCor_body <- function(processed,
                  ranPars=NULL, ## all dispersion and correlation params ideally provided through ranPars
                  distMatrix,uniqueGeo=NULL,adjMatrix,corrMatrix,
                  verbose=c(warn=TRUE,trace=FALSE,summary=FALSE),control.dist=list(),
                  ...) { 
  dotlist <- list(...)
  ################# 
  data <- processed$data
  predictor <- processed$predictor 
  spatial.terms <- findSpatial(predictor)
  spatial.model <- spatial.terms[[1L]] 
  if ( ! is.null(spatial.model)) {
    corr.model <- as.character(spatial.model[[1L]]) 
  } else {
    if ( ! missing(corrMatrix)) {
      mess <- pastefrom("corrMatrix argument despite no corrMatrix term in formula:",prefix="(!) From ")
      message(mess)
      stop("This syntax is obsolete; add a corrMatrix(...) term in the formula.")
    } ## ELSE more generic message: 
    stop("Call to 'HLCor' without a spatial term in the formula is suspect.")
  }
  ## convert back ranPars to canonical scale:
  rpblob <- canonizeRanPars(ranPars=ranPars,corr.model=corr.model) 
  ranPars <- rpblob$ranPars
  trueCorrpars <- rpblob$trueCorrpars
  rho <- ranPars$rho
  #
  coordinates <- NULL
  test.in <- FALSE
  ### ensure LMatrix in predictor: 
  ## if it is currently absent, first provide corr matrix or its symSVD, from which Lunique will be computed using designL.from.Corr
  if (is.null(Lunique <- attr(predictor,"LMatrix"))) { 
    symSVD <- NULL
    if (corr.model %in% c("adjacency","ar1")) { ## "ar1" != "AR1" is a tempo name for a futur generic model  
      if ( missing(adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
      ## no nugget in the adjacency model... ## (use additional ranef instead)      
      symSVD <- attr(adjMatrix,"symSVD")
      if (is.null(symSVD) && identical(attr(ranPars,"type")$rho,"var")) { ## can occur in direct call of HLCor ## identical() handles NULL args
        if (isSymmetric(adjMatrix)) {
          symSVD <- selfAdjointWrapper(adjMatrix)
          attr(adjMatrix,"symSVD") <- symSVD
        }             
      }
      if (is.null(symSVD)) {
        corrm <- solve(diag(nrow(adjMatrix))-rho*(adjMatrix))
      } else {
        symSVD$adjd <- symSVD$d
        symSVD$d <- 1/(1-rho*symSVD$d) ## from adjMatrix to correlation matrix
      }
    }  else if (corr.model %in% c("SAR_WWt")) { ## "ar1" != "AR1" is a tempo name for a futur generic model  
      if ( missing(adjMatrix) ) stop("missing 'adjMatrix' for adjacency model")
      UDU. <- attr(adjMatrix,"UDU.")
      if (is.null(UDU.)) {
        corrm <- solve(diag(nrow(adjMatrix))-rho*(adjMatrix))
      } else {
        corrm <- UDU.$u %*% sweep(UDU.$u.,MARGIN=1,1/(1-rho*UDU.$d),`*`) 
      }
      corrm <- tcrossprodCpp(corrm)
    }  else if (corr.model=="AR1") {
      coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(data))
      uniqueGeo <- calcUniqueGeo(data=data[,coordinates,drop=FALSE])
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        stop("HLCor code should be allowed again to handle blockDiag objects")
        #scaled.dist <- as.blockDiag.bar(spatial.model[[2]],formula,data=uniqueGeo)
        #test.in <- TRUE
      } else scaled.dist <- proxy::dist(uniqueGeo)
      corrm <- trueCorrpars$ARphi^scaled.dist
    } else  if (corr.model %in% c("Matern")) {
      txt <- paste(spatial.model[[2]][[3]]) ## the RHS of the ( . | . ) 
      if (length(grep("%in%",txt))>0) {
        stop("(!) Matern( . | <coord> %in% <grp>) is not yet handled.")
        test.in <- TRUE ## should be useful when this case will be handled
      } 
      ## in a typical call from corrHLfit the following test should be FALSE because uniqueGeo and maybe distMatrix should have been precomputed
      if ((length(rho)>1 || missing(distMatrix)) && is.null(uniqueGeo)) { ## all cases where we need uniqueGeo
        coordinates <- extract.check.coords(spatial.model=spatial.model,datanames=names(data))
        uniqueGeo <- calcUniqueGeo(data=data[,coordinates,drop=FALSE]) ## keeps the names of first instances of the coordinates in data
      } 
      ## then compute scaled distances from unscaled info, for HLfit call
      msd.arglist <- list(rho = rho)
      msd.arglist$`dist.method` <- control.dist$`dist.method` ## may be NULL
      if (length(rho)>1L) {
        msd.arglist <- c(msd.arglist,list(uniqueGeo=uniqueGeo))
        msd.arglist$`rho.mapping` <- control.dist$`rho.mapping` ## may be NULL
      } else {
        if ( missing(distMatrix)) { 
          dist.arglist <- list(x=uniqueGeo)
          dist.arglist$method <- control.dist$dist.method ## may be NULL
          distMatrix <- do.call(proxy::dist,dist.arglist)
        }
        msd.arglist <- c(msd.arglist,list(distMatrix=distMatrix))
      }
      corrm <- do.call("make_scaled_dist",msd.arglist)
      ## at this point is a single location, corrm should be dist(0) and make_scaled_dist was modified to that effect
      if ( nrow(corrm)>1 ) { ## >1 locations
        norho <- trueCorrpars; norho$rho <- NULL ## because the MaternCorr input will be an already scaled distance 'corrm'
        corrm <- do.call(MaternCorr,args=c(norho,list(corrm)))        
      } 
    } else if (corr.model== "corrMatrix") {
      if (missing(corrMatrix)) {
        mess <- pastefrom("missing(corrMatrix) argument despite corrMatrix term in formula.",prefix="(!) From ")
        stop(mess)
      } ## ELSE:
      corrm <- corrMatrix
      Lunique <- attr(corrMatrix,"LMatrix") ## will typically be NULL, but super-users ;-) may have provided it
    } 
    if (verbose["trace"] && length(trueCorrpars)>0) print(unlist(trueCorrpars))
    ## call designL.from.Corr if Lunique not available
    if (is.null(Lunique)) { ## test FR 11/2013 ## modif 2015/04. Noter un calcul de Lunique ci dessus
      if ( ! is.null(symSVD)) {
        Lunique <- try(designL.from.Corr(symSVD=symSVD))
      } else { ## corrm must exist
        argsfordesignL <- dotlist[intersect(names(dotlist),names(formals(designL.from.Corr)))] 
        if (processed$HL[1L]=="SEM") argsfordesignL$try.chol <- FALSE
        if (inherits(corrm,"dist")) {
          corrm <- as.matrix(corrm)
          diag(corrm) <- 1L ## always a correlation matrix
        }
        Lunique <- try(do.call(designL.from.Corr,c(list(m=corrm),argsfordesignL)))
      }
      if (inherits(Lunique,"try-error")) { 
        print("correlation parameters were:") ## makes sense if designL.from.Corr already issued some warning
        print(unlist(trueCorrpars))    
        stop()
      }
    }
    attr(predictor,"%in%") <- test.in
    attr(Lunique,"corr.model") <- corr.model
    attr(Lunique,"ranefs") <- unlist(lapply(spatial.terms,DEPARSE)) ## essentiel pour la construction de ZAL!
    if ( corr.model=="adjacency"
         && ! is.null(attr(ranPars,"type")) ## ie through corrHLfit call
         && "var" %in% attr(ranPars,"type")$rho ## then not a call for fixed rho => estim of rho within HLfit through SEM or augm GLM
         ) { ## then define ZA.L as ZA. U(adjacency matrix)
      Lunique[] <- attr(Lunique,"symsvd")$u ## "[]" keeps attributes
    }
    attr(predictor,"LMatrix") <- Lunique
  }
  processed$predictor <- predictor
  ###
  HLFormals <- names(formals(HLfit))
  good_dotnames <- intersect(names(dotlist),HLFormals)
  if (length(good_dotnames)>0L) {
    HL.info <- dotlist[good_dotnames]
  } else HL.info <- list()
  ## all printing in HLfit is suppressed by default
  HL.info$verbose <- verbose #[intersect(names(verbose),c("warn","trace","summary","SEM"))] 
  HL.info$processed <- processed
  ## convert ranPars to ranFix + init.HLfit
  ## allows log and not log:
  varNames <- names(which(attr(ranPars,"type")=="var"))
  HL.info$init.HLfit[varNames] <- ranPars[varNames] ## inherits values from corrHLfit(...,init.HLfit(...))
  fixNames <- setdiff(names(ranPars),varNames) 
  if (!is.null(fixNames)) { ## could be NULL for corrMatrix case
    ranFix <- ranPars[fixNames] ## 11/2014 as there is no other source for ranFix
    typelist <- list() 
    typelist[fixNames] <- "fix" 
    if (!is.null(rPtype <- attr(ranPars,"type"))) { ## it may not exist, or elements may be "fix" or "outer"
      typelist[names(rPtype)] <- rPtype
    }
    attr(ranFix,"type") <- typelist 
    HL.info$ranFix <- ranFix
  }
  hlfit <- do.call("HLfit",HL.info) ## with a _list_ of arguments -> do.call ## perhaps should create a list of unevaluated arguments ???? 
  if ( ! is.null(hlfit$error)) {
    errfile <- generateFileName("HLfitCall")
    errfile <- paste(errfile,".RData",sep="")
    save(HL.info,file=errfile)
    mess <- pastefrom("'do.call(HLfit,HL.info)' failed:",prefix="(!) From ")
    message(mess)
    message(hlfit$error)
    message("'HL.info' is saved in the ",errfile," file",sep="")
    stop("I exit.")
  } ## ELSE:
  hlfit$control.dist <- control.dist
  attr(hlfit,"info.uniqueGeo") <- uniqueGeo ## more spatial info is in the hlfit$predictor's attributes (Lunique = corrmat^1/2, and ZALMatrix)
  if (corr.model %in% c("Matern")) attr(hlfit,"msd.arglist") <- msd.arglist ## more organized, easier to reuse. 
  # particulary for $rho.mapping:
  #  should be NULL if length rho = 1 or if original control.dist$rho.mapping was NULL
  ## FR->FR but info.uniqueGeo more general (eg AR1) -> a revoir
  hlfit$call <- "$call removed by HLCor. Consider the 'HLCorcall' attribute instead." ## instead of the $call with evaluated arguments
  return(hlfit) ## 
}


## wrapper for HLCor, suitable input and output for optimization
`HLCor.obj` <- function(ranefParsVec,skeleton,HLCor.obj.value="p_bv",traceFileName=NULL,...) { ## name of first arg MUST differ from names in dotlist...
  mc <- match.call(expand.dots=TRUE) ## (1) expand.dots added 11/04/2014 for the multinomial... eval 
  if (is.null(processed <- mc$processed)) {
    stop("Call to HLCor.obj() without a 'processed' argument is invalid")
  } else { ## 'processed' is available
    multiple <- attr(processed,"multiple")
    if ( ( ! is.null(multiple)) && multiple)  { ## "multiple" processed list 
      ## RUN THIS LOOP and return
      fitlist <- lapply(seq_len(length(processed)),function(it){
        locmc <- mc
        locmc[[1L]] <- as.name("HLCor.obj") ## replaces "f" !
        locmc$ranefParsVec <- ranefParsVec ## replaces "arg" !
        locmc$processed <- processed[[it]] ## The data are in processed !
        locmc$distMatrix <- mc$distMatrix[[it]] ## but the matrices are not HLfit args hence not in processed ! 
        locmc$uniqueGeo <- mc$uniqueGeo[[it]]
        eval(locmc)
      }) ## a pure list of HLCor objects
      resu <- sum(unlist(fitlist))
      if (is.character(traceFileName)) {
        verif <- paste("#global:",ranefParsVec,resu) 
        write(verif,file=traceFileName,append=T) ## the file is unlink'ed in corrHLfit()  
      }
      return(resu)
    } else { ## there is one processed for a single data set 
      family <- processed$family
      data <- processed$data
    }
  }
  
  HLCor.formals <- names(formals(HLCor))
  names_formals_HLfit <- names(formals(HLfit))
  designL.formals <- names(formals(designL.from.Corr))
  makescaled.formals <- names(formals(make_scaled_dist))
  HLnames <- (c(HLCor.formals,names_formals_HLfit,designL.formals,makescaled.formals))  ## cf parallel code in corrHLfit
  HLCor.call <- mc[c(1,which(names(mc) %in% HLnames))] ## keep the call structure
  HLCor.call[[1L]] <- quote(spaMM::HLCor)
  forGiven <- relist(ranefParsVec,skeleton) ## given values of the optimized variables 
  ## ... relist keeps the RHOMAX... attributes from the skeleton, but the partial copy into ranPars does not. 
  HLCor.call$ranPars[names(forGiven)] <- forGiven ## do not wipe out other fixed, non optimized variables
  attr(HLCor.call$ranPars,"RHOMAX") <- attr(skeleton,"RHOMAX")
  attr(HLCor.call$ranPars,"NUMAX") <- attr(skeleton,"NUMAX")
  types <- attr(skeleton,"type")
  attr(HLCor.call$ranPars,"type")[names(types)] <- types
  if (is.character(traceFileName)) {
    if(.spaMM.data$options$TRACE.UNLINK) unlink("HLCor.call*.RData")
    zut <- paste(ranefParsVec,collapse="")  
    save(HLCor.call,file=paste("HLCor.call",zut,".RData",sep="")) ## for replicating the problem
  }
  hlfit <- eval(HLCor.call)
  aphls <- hlfit$APHLs
  resu <- aphls[[HLCor.obj.value]]
  if (is.character(traceFileName)) {
    readable <- unlist(canonizeRanPars(ranPars=forGiven,corr.model=mc$`corr.model`,checkComplete=FALSE)$ranPars) 
    verif <- c(unlist(aphls),hlfit$lambda,hlfit$phi,readable,ranefParsVec) ## hlfit$phi may be NULL
    write(verif,file=traceFileName,ncolumns=length(verif),append=TRUE) ## the file is unlink'ed in corrHLfit()  
  }
  return(resu) #
}


