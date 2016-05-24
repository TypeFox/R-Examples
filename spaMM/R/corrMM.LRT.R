`Bartlett.robust` <- function(LRTobject,robust=TRUE,verbose=F) {
      bootLRTS <- with(LRTobject,2*(bootreps[,1]-bootreps[,2]))## full -null but can be p_v or p_bv
      # plot(qchisq(ppoints(zut),1),sort(zut)) ## QQplot, MASS p. 108
      filter <- bootLRTS[bootLRTS>-1e-08]
      resu <- list()
      if (robust) {
        ## finds the df of the distribution by robust regression (MM: MASS p. 161) of QQplot, 
        robustMean <- rlm(sort(filter)~qchisq(ppoints(filter),1)-1,method="MM",maxit=200)$coefficients[[1]] 
        ## [[1]] to remove name which otherwise finishes as a rowname in a subsequent dataframe       
        # plot(qchisq(ppoints(filter),1),sort(filter)) ## QQplot, MASS p. 108
        # points(qchisq(ppoints(filter),1),robustMean * qchisq(ppoints(filter),1),pch=".")   
        if (inherits(robustMean,"try-error")) {
          resu$robustMean <- NA
          mess <- pastefrom("problem in computation of robustMean.")
          warning(mess)
        } else resu$robustMean <- robustMean
      }
      resu$meanPosbootLRT <- mean(filter)
      resu$nPosbootLRT <- length(filter) 
      if (verbose) print(unlist(resu))
      return(resu)
}

## fixedLRT is a safe interface for performing tests as described in Ecography paper.
## it does not allow the profiling procedure in corrMM.LRT
fixedLRT <- function(null.formula,formula,data,HLmethod,REMLformula=NULL,boot.repl=0,
                        control=list(),control.boot=list(),...) {  ## since corrMM.LRT is not doc'ed, REMLformula=NULL,boot.repl=0 cannot go into '...' 
  if (missing(null.formula)) stop("'null.formula' argument is missing, with no default.")
  if (missing(formula)) stop("'formula' argument is missing, with no default.")
  if (missing(data)) stop("'data' argument is missing, with no default.")
  if (missing(HLmethod)) stop("'HLmethod' argument is missing, with no default.")
  if (! is.null(control$profiles)) {
    stop("'fixedLRT' does not allow 'control$profiles'.")
  }
  ## see 'lm' code for template
  mc <- match.call(expand.dots = TRUE)
  ## other possible settings, through iterative fits
  #  mc$init.corrHLfit$lambda <- NULL
  #  mc$init.corrHLfit$phi <- NULL
  ## we have a potential backward compatiblity problem, since the simulation scripts 
  ## for the Ecography paper assume that the package automatically interpret the model as spatial, even if findSpatial returns NULL
  ## and we no longer want such a behaviour
  ## but fixedLRT is not used in these scripts, so it can make a different assumption
  spatial <- findSpatial(formula)
  if ( ! is.null(spatial)) {
    mc$method <- "corrHLfit" ## useful for spaMMLRT call but not for corrMM.LRT where it is the default
    ## both will use p_v for the optim steps, we need to distinguish whether some REML correction is used in iterative algo :
    if ( HLmethod %in% c("ML","PQL/L","SEM") || substr(HLmethod,0,2) == "ML") {
      mc[[1L]] <- as.name("spaMMLRT") ## does not (yet) handles well other HLmethod's  when eg init.corrHLfit contains lambda
      ## there's no profile etc in spaMMLRT... 
    } else { ## EQL, REPQL or REML variants: profiles then not allowed within corrMM.LRT!
      mc[[1L]] <- as.name("corrMM.LRT") ## corrMM.LRT methods and its options below are best frozen to their v1.0 state
      mc$control<-list(profiles=0,prefits=FALSE) ## default values in call by fixedLRT. corrMM.LRT further has default restarts=TRUE and maxit=1
      mc$control[names(control)]<-control ## overrides with user values
      mc$control.boot <- control.boot ## default values in call by fixedLRT are those of corrMM.LRT ie prefits=FALSE,profiles=0. We can directly copy user values. 
    }
  } else {
    mc[[1L]] <- as.name("spaMMLRT") 
    ## No profiles, maxit, restarts, prefits
    if (is.null(mc$corrMatrix)) { ## neither explicit spatial nor corrMatrix -> HLfit
        mc$method <- "HLfit"
      } else {
      ## corrMatrix -> we need to use HLCor
        mc$method <- "HLCor"
      }                 
  }  
  eval(mc, parent.frame())
}


spaMMLRT <- function(null.formula=NULL,formula,
                     null.disp=list(),REMLformula=NULL,
                     method="corrHLfit",boot.repl=0,
                     ## currently trace always false; this is not an argument t be forwarded as is to corrHLfit! 
                     trace=FALSE, ## T means lead to calls of corrHLfit(... trace=list(<file name>,<over/append>))
                     verbose=c(trace=FALSE,warn=NA,summary=FALSE),  
                     ...) {
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- FALSE ## will be unconditionally ignored by the final fit in corrHLfit  
  if (is.na(verbose["summary"])) verbose["summary"] <- FALSE ## this is for HLCor
  dotlist <-list(...)
  ## birth pangs :
  if ("predictor" %in% names(dotlist)) {
    stop("'spaMMLRT' called with 'predictor' argument which should be 'formula'" )
  }
  if ("null.predictor" %in% names(dotlist)) {
    stop("'spaMMLRT' called with 'null.predictor' argument which should be 'null.formula'" )
  }
  ## here we makes sure that *predictor variables* are available for all data to be used under both models
  data <- dotlist$data
  if ( inherits(data,"list")) {
    data <- lapply(data,function(dt) {
      null.validdata <- validData(formula=null.formula[-2],resid.formula=dotlist$resid.formula,data=dt) ## will remove rows with NA's in required variables
      full.validdata <- validData(formula=formula[-2],resid.formula=dotlist$resid.formula,data=dt) ## will remove rows with NA's in required variables
      dt[intersect(rownames(null.validdata),rownames(full.validdata)),,drop=FALSE]     
    })
  } else {
    null.validdata <- validData(formula=null.formula[-2],resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    full.validdata <- validData(formula=formula[-2],resid.formula=dotlist$resid.formula,data=data) ## will remove rows with NA's in required variables
    data <- data[intersect(rownames(null.validdata),rownames(full.validdata)),,drop=FALSE]     
  }  
  dotlist$data <- data
  predictor <- formula   
  if (! inherits(formula,"predictor")) predictor <- Predictor(formula)
  null.predictor <- null.formula   
  if (! inherits(null.formula,"predictor")) null.predictor <- Predictor(null.formula)
  form <- predictor
  if (!is.null(dotlist$LamFix)) {
    dotlist$ranFix$lambda <- dotlist$LamFix
    dotlist$LamFix <- NULL
  }
  if (!is.null(dotlist$PhiFix)) {
    dotlist$ranFix$phi <- dotlist$PhiFix
    dotlist$PhiFix <- NULL
  }  
  dotlist$formula <- predictor 
  dotlist$verbose <- verbose 
  if (method=="corrHLfit") {
    if (is.null(dotlist$objective)) { ## implements default objective function
      if ( ! is.null(null.predictor)) {dotlist$objective <- "p_v"} else {dotlist$objective <- "p_bv"}
    }
  } else dotlist$objective <- NULL
  ## do not change dotlist afterwards !
  fullm.list <- dotlist
  nullm.list <- dotlist
  fullm.list$formula <- predictor
  #### "limited use, for initializing bootstrap replicates:"
  which.iterative.fit <-character(0)
  which.optim.fit <-character(0)
  if ( ! is.null(dotlist$init.corrHLfit$lambda) ) {
    which.optim.fit <- c(which.optim.fit,"lambda")
  } else if ( is.null (fullm.list$ranFix$lambda)) which.iterative.fit <- c(which.iterative.fit,"lambda")
  if (is.null(fullm.list$family) ## (=gaussian)
      || fullm.list$family$family %in% c("gaussian","Gamma")) {
    if ( ! is.null(dotlist$init.corrHLfit$phi) ) { ## if user estimates it by ML...
      which.optim.fit <- c(which.optim.fit,"phi")
    } else which.iterative.fit <- c(which.iterative.fit,"phi")
  }
  if (method=="corrHLfit") {
    if ( ! "rho" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"rho")
    if ( ! "nu" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"nu")
    if ( ! "ARphi" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"ARphi") # mpf
    if ( ! "Nugget" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"Nugget") ## by default not operative given later code and init.optim$Nugget is NULL
  }
  if ( ! is.null(null.predictor)) { ## ie if test effet fixe
    testFix <- T
    if (dotlist$HLmethod =="SEM") {
      test.obj <- "logLapp"
    } else test.obj <- "p_v"
    ## check fullm.list$REMLformula, which will be copied into nullm in all cases of fixed LRTs
    if (dotlist$HLmethod %in% c("ML","PQL/L","SEM") || substr(dotlist$HLmethod,0,2) == "ML") {
      fullm.list$REMLformula <- NULL
      nullm.list$REMLformula <- NULL
    } else { ## an REML variant
      if (is.null(REMLformula)) { ## default
        fullm.list$REMLformula <- NULL ## HLfit will reconstruct proper REML from this
        nullm.list$REMLformula <- fullm.list$formula ## 
      } else { ## allows alternative choice, but still the same *both* the full and the null fit
        fullm.list$REMLformula <- REMLformula
        nullm.list$REMLformula <- REMLformula
      }
      namesinit <- names(fullm.list$init.corrHLfit)
      namesinit <- setdiff(namesinit,c("rho","nu","Nugget","ARphi"))
      len <- length(namesinit)
      if ( len>0) {
        if (len > 1) namesinit <- paste(c(paste(namesinit[-len],collapse=", "),namesinit[len]),collapse=" and ")
        message("Argument 'init.corrHLfit' is used in such a way that")
        message(paste("  ",namesinit," will be estimated by maximization of p_v.",sep=""))
        message("  'REMLformula' will be inoperative if all dispersion")
        message("  and correlation parameters are estimated in this way.")
      }
    }
    nullm.list$formula <- null.predictor
  } else if ( length(null.disp)>0 ) { ## test disp/corr param
    testFix <- F
    test.obj <- "p_bv"
    #
    mess <- pastefrom("changes to objective fn required to renullfit or refull fit computation.")
    stop(mess)
    #
    namescheck <- names(fullm.list$lower)
    namescheck <- namescheck[ ! namescheck %in% names(null.disp)]  
    nullm.list$lower <- nullm.list$lower[namescheck]
    nullm.list$upper <- nullm.list$upper[namescheck]
    nullm.list$ranFix <- c(nullm.list$ranFix,null.disp) ## adds fixed values to any preexisting one
  } else testFix <- NA
  
#  trace.info <- NULL
  fullfit <- do.call(method,fullm.list)
  nullfit <- do.call(method,nullm.list)
  if (testFix) {df <- length(fullfit$fixef)-length(nullfit$fixef)} else {df <- length(null.disp)}
  if (df<0) {
    tmp <- fullfit
    fullfit <- nullfit
    nullfit <- tmp
  }
  if (inherits(fullfit,"HLfitlist")) {
    fullL <- attr(fullfit,"APHLs")[[test.obj]]
  } else fullL <- fullfit$APHLs[[test.obj]]
  if (inherits(nullfit,"HLfitlist")) {
    nullL <- attr(nullfit,"APHLs")[[test.obj]]
  } else nullL <- nullfit$APHLs[[test.obj]]
  LRTori <- 2*(fullL-nullL)
  ## BOOTSTRAP
  if ( ! is.na(testFix)) {
    if (boot.repl>0) {
      bootlist <- dotlist ## copies (full)formula and (optionally) ranFix
      bootlist <- c(bootlist,list(null.formula=null.predictor,null.disp=null.disp,REMLformula=REMLformula,method=method)) ## unchanged user REMLformula forwarded
      bootlist$verbose <- c(trace=FALSE,summary=FALSE)
      bootlist$trace <- FALSE 
      bootlist$boot.repl <- 0 ## avoids recursive call of bootstrap
      all.estim.ranvars <- c(which.optim.fit) ## FR->FR to be modified if optimFits are reintroduced, see corresponding code in corrMM.LRT 
      for (st in all.estim.ranvars) {
        if (dotlist$HLmethod=="corrHLfit") {
          if ( st %in% names(nullfit$corrPars)) {
            bootlist$init.corrHLfit[st] <- nullfit$corrPars[st]
          } else if ( st %in% names(nullfit)) {
            bootlist$init.corrHLfit[st] <- nullfit[st] ## handled in dotlist by corrHLfit, cf notes 090113
          } ## it's also possible that st is nowhere (10/2013: currently for Nugget)
        } else {
          if ( st %in% names(nullfit)) {
            bootlist$init.HLfit[st] <- nullfit[st] ## 
          } 
        }
      }        
      bootreps<-matrix(,nrow=boot.repl,ncol=2) 
      colnames(bootreps) <- paste(c("full.","null."),test.obj,sep="")
      msg <- "bootstrap replicates: "
      msglength <- nchar(msg)
      cat(msg)
      t0 <- proc.time()["user.self"]
      simbData <- nullfit$data
      if (tolower(nullfit$family$family)=="binomial") {
        nform <- attr(nullfit$predictor,"oriFormula")  
        if (is.null(nform)) stop("a 'predictor' object must have an 'oriFormula' member.")
        if (paste(nform[[2L]])[[1L]]=="cbind") {
          ## We have different possible (exprL,exprR) arguments in cbind(exprL,exprR), 
          ## but in all case the predictor is that of exprL and exprR is $weights- exprL. We standardize: 
          nposname <- makenewname("npos",names(data))
          nnegname <- makenewname("nneg",names(data))
          nform <- paste(nform)
          nform[2L] <- paste("cbind(",nposname,",",nnegname,")",sep="")
          bootlist$null.formula <- as.formula(paste(nform[c(2,1,3)],collapse=""))
          fform <- paste(attr(fullfit$predictor,"oriFormula"))
          fform[2L] <- nform[2L]
          bootlist$formula <- as.formula(paste(fform[c(2,1,3)],collapse=""))
          cbindTest <- TRUE
        } else cbindTest <- FALSE
      } else cbindTest <- FALSE
      ## the data contain any original variable not further used; e.g original random effect values in the simulation tests  
      thisFnName <- as.character(sys.call()[[1]]) ## prevents a bug when we change "this" function name
      for (ii in 1:boot.repl) {
        locitError <- 0
        repeat { ## for each ii!
          newy <- simulate(nullfit) ## cannot simulate all samples in one block since some may not be analyzable  
          if (cbindTest) {
            simbData[[nposname]] <- newy
            simbData[[nnegname]] <- nullfit$weights - newy
          } else {simbData[[as.character(nullfit$predictor[[2L]])]] <- newy} ## allows y~x syntax for binary response
          bootlist$data <- simbData
          bootrepl <- try(do.call(thisFnName,bootlist)) ###################### CALL ##################
          if (! inherits(bootrepl,"try-error") ) { ## eg separation in binomial models... alternatively, test it here (require full and null X.pv... )
            if (inherits(bootrepl$fullfit,"HLfitlist")) {
              fullL <- attr(bootrepl$fullfit,"APHLs")[[test.obj]]
            } else fullL <- bootrepl$fullfit$APHLs[[test.obj]]
            if (inherits(bootrepl$nullfit,"HLfitlist")) {
              fullL <- attr(bootrepl$nullfit,"APHLs")[[test.obj]]
            } else fullL <- bootrepl$nullfit$APHLs[[test.obj]]
            bootreps[ii,] <- c(bootrepl$fullfit$APHLs[[test.obj]],bootrepl$nullfit$APHLs[[test.obj]])
            break ## replicate performed, breaks the repeat
          } else { ## there was one error
            locitError <- locitError + 1
            if (locitError>10) { ## to avoid an infinite loop
              stop("Analysis of bootstrap samples fails repeatedly. Maybe no statistical information in them ?")
            } ## otherwise repeat!
          }
        } 
        tused <- proc.time()["user.self"]-t0
        ttotal <- tused* boot.repl/ii
        if (interactive()) {
          for (bidon in 1:msglength) cat("\b")
          msg <- paste("Estimated time remaining for bootstrap: ",signif(ttotal-tused,2)," s.",sep="")
          msglength <- nchar(msg)
          cat(msg)
        } else {
          cat(ii);cat(" ")
          if ((ii %% 40)==0L) cat("\n")
        }
      }
      cat("\n")
    } ## end main bootstrap loop
  } else { ## nothing operativ yet
    bootreps<-matrix(,nrow=boot.repl,ncol=length(unlist(fullfit$APHLs))) 
    colnames(bootreps) <- names(unlist(fullfit$APHLs))
    ## more needed here ?
  }
  ## prepare output
  if ( ! is.na(testFix)) {
    if (testFix) {df <- length(fullfit$fixef)-length(nullfit$fixef)} else {df <- length(null.disp)}
    resu <- list(fullfit=fullfit,nullfit=nullfit)
    LRTinfo <- list(df=df,LRTori = LRTori)
    basicLRT <- data.frame(LR2=LRTori,df=df,pvalue=1-pchisq(LRTori,df=df))
    if (boot.repl>0) {
      bootdL <- bootreps[,1]-bootreps[,2]
      meanbootLRT <- 2*mean(bootdL)  
      LRTinfo$rawPvalue <- (1+sum(bootdL>=LRTori/2))/(boot.repl+1) ## DavisonH, p.141
      LRTinfo$meanbootLRT <- meanbootLRT
      LRTinfo$bootreps <- bootreps
      LRTinfo$LRTcorr <- LRTori*df/meanbootLRT
    }
  } else {
    resu <- list(fullfit=fullfit)
    LRTinfo <- list()
    basicLRT <- list()
  }
  #  LRTinfo$trace.info <- trace.info 
  ##  resu$LRTinfo <- LRTinfo ## pas compatible avec hglmjob.R...
  resu <- c(resu,LRTinfo) ## loses the sublist structure, which wouldnot be compatible with hglmjob.R...  
  ## this keeps the dataframe structure of basicLRT:
  resu <- c(resu,list(basicLRT = basicLRT)) ## as in anova and LRT and their summaries  
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}


`corrMM.LRT` <- function(null.formula=NULL,formula,
                       null.predictor=null.formula,predictor=formula, ## simple back compat code...
                       null.disp=list(),REMLformula=NULL,
                       method="corrHLfit",boot.repl=0,
                       which.iterative=c(),
                       trace=FALSE, ## T means lead to calls of corrHLfit(... trace=list(<file name>,<over/append>))
                       control=list(), ## profiles=Inf,prefits=T,optimFits,restarts,maxit...
                       control.boot=list(), ## prefits=F,optimFits,profiles=0
                       verbose=c(trace=FALSE,warn=NA,summary=FALSE),  
                       test.obj=NULL, ## 01/2014 to override default behaviour (experimental)  
                       ...) {
  if (is.na(verbose["trace"])) verbose["trace"] <- FALSE
  if (is.na(verbose["warn"])) verbose["warn"] <- FALSE ## will be unconditionally ignored by the final fit in corrHLfit  
  if (is.na(verbose["summary"])) verbose["summary"] <- FALSE ## 
  ##################################################                     
  ## methods handled:
  ## RE<...> (ie. 'REML', 'RE(0,1)'...) => REML fits using the iterative algorithm, both with design matrix from full models; a bit slow...
  ## Other methods: no REML correction but the same REML fits are still possible to provide starting values for other fits
  ## if an init.corrHLfit argument is provided, then there is no REML estimation of the given parameter, even if the given method is RE<...>
  ## 12/2013: PQL now identified to REPQL. At some stage there was 'monPQL' which was ML(0,1)
  ## my HL(...) were ML(...) because of init.corrHLfit, but would otherwise be RE(...)
  ## FR->FR cmmentaire non daté << v1.1 suggerait problèmes persistants ML(...) + REML formula. Pt etre pas un cas int
  ##################################################                     
  profiles <- control$profiles
  if (is.null(profiles)) profiles <- Inf
  restarts <- control$restarts
  if (is.null(restarts)) restarts <- TRUE
  maxit <- control$maxit
  if (is.null(maxit)) maxit <- 1
  REMLfits <- control$REMLfits
  if ( ! is.null(REMLfits)) { ## back compat code for change of name
    control$prefits <- REMLfits
    control$REMLfits <- NULL
  }
  # added 10/2013 probably in anticipation for direct uses of corrMM.LRT
  # the prefits are F by default
  # I should be able to use boot.control directly 
  bootcontrol <- list(prefits=FALSE,profiles=0) ## corrMM/LRT defaults (optimFits defaults controlled elsewhere)
  bootcontrol[names(control.boot)] <- control.boot ## overrriden by user input
  #
  bootFix <- control$bootFix
  if (is.null(bootFix)) bootFix <- c()
  if (! inherits(predictor,"predictor")) predictor <- Predictor(predictor)
  form <- predictor
  dotlist <-list(...)
  if (!is.null(dotlist$LamFix)) {
      dotlist$ranFix$lambda <- dotlist$LamFix
    dotlist$LamFix <- NULL
  }
  if (!is.null(dotlist$PhiFix)) {
    dotlist$ranFix$phi <- dotlist$PhiFix
    dotlist$PhiFix <- NULL
  }  
  dotlist$formula <- predictor 
  dotlist$verbose <- verbose 
  if (method=="corrHLfit") {
    if (is.null(dotlist$objective)) { ## implements default objective function
      if ( ! is.null(null.predictor)) {dotlist$objective <- "p_v"} else {dotlist$objective <- "p_bv"}
    }
  } else dotlist$objective <- NULL
  ## do not change dotlist afterwards !
  trace.info <- NULL
  nullm.list <- dotlist
  fullm.list <- dotlist
  ####
  which.iterative.fit <-character(0)
  which.optim.fit <-character(0)
  if ( ! is.null(dotlist$init.corrHLfit$lambda) ) { ## if user estimates it by ML...
    which.optim.fit <- c(which.optim.fit,"lambda")
  } else if ( is.null (fullm.list$ranFix$lambda)) which.iterative.fit <- c(which.iterative.fit,"lambda")
  if (is.null(fullm.list$family) ## (=gaussian)
      || fullm.list$family$family %in% c("gaussian","Gamma")) {
    if ( ! is.null(dotlist$init.corrHLfit$phi) ) { ## if user estimates it by ML...
      which.optim.fit <- c(which.optim.fit,"phi")
    } else which.iterative.fit <- c(which.iterative.fit,"phi")
  }
  if ( ! "rho" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"rho")
  if ( ! "nu" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"nu")
  if ( ! "Nugget" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"Nugget") ## by default not operative
  if ( ! "ARphi" %in% names(dotlist$ranFix)) which.optim.fit <- c(which.optim.fit,"ARphi") ## hmf
######## Determines iterativeFits, optimFits and the corresponding REML formulas
  ## *by default* the iterative fits are more REML than the optim fits; true REML is possible only through the iterative fits.
  #### optimREMLformula is a (false, typically ML) REML formula for nullfit, fullfit, not for REML fits
  if (is.null(dotlist$HLmethod)) stop("Argument 'HLmethod' missing with no default in corrMM.LRT") ## Devrait pas venir de la dotlist... codage bouseux mais bon.
  if (dotlist$HLmethod %in% c("ML","PQL/L") || substr(dotlist$HLmethod,0,2) == "ML") { ## 
    iterativeFits <- control$prefits ## performs an iterative prefit for dispersion params that have afterwards fitted using optim
    optimFits <- control$optimFits ## 
    if (is.null(optimFits)) optimFits <- TRUE 
    if (optimFits) optimREMLformula <- NULL  ## for ML preprocess will reject any non null input REML formula and reconstruct it internally  
  } else { ## RE, HL, EQL... will skip the loop; for consistency with RE in v1.0, and by lack of interest otherwise.  
    iterativeFits <- TRUE
    optimFits <- FALSE ## hence profiles not possible
  } 
  if ( (!optimFits) && profiles>0 ) {stop("Profiles not possible with HLmethod other than 'ML...' or 'PQL/L'.")}
  ## v1.0 had:
#  if (dotlist$HLmethod %in% c("ML","PQL/L") || substr(dotlist$HLmethod,0,2) == "ML") { ## 
#    iterativeFits <- control$prefits 
#    optimFits <- control$optimFits ## 
#    if (is.null(optimFits)) optimFits <- TRUE 
#    if (optimFits) optimREMLformula <- NULL  ## for ML preprocess will reject any non null input REMLformula and reconstruct it internally  
#  } else if (substr(dotlist$HLmethod,0,2) %in% c("RE")) { ## explicit request for REML; only iterativeFits should be run and optimREMLformula shoul not be used 
#    iterativeFits <- TRUE
#    optimFits <- FALSE
#  } else {  ## HL<...>: or EQL. In corrMM.LRT, this was interpreted as ML fits. But this is reversed as non-standard in v1.1.   
#    iterativeFits <- control$prefits 
#    bars <- findbarsMM(form)  ## extract random effects
#    lhs <- paste(form[[2]]) ## extract response
#    ## build formula with only random effects
#    optimFits <- control$optimFits ## 
#    if (is.null(optimFits)) optimFits <- TRUE 
#    if (optimFits) optimREMLformula <- as.formula(paste(lhs,"~",paste(paste("(",as.character(bars),")"),collapse="+"))) ## ML in iterative fit
#  }
  ####
  if (is.null(iterativeFits)) { ## ie default for non RE<...>
    if (length(which.iterative)==0L) {
      iterativeFits <- TRUE ## we can construct a non trivial iterative fit that usefully can be used to initiate the optim fits
    } else {
      iterativeFits <- FALSE ## then the optim fits are iterative wrt to lambda in *G*LMMs; then we usually don't need the REML fits to find a good starting lambda, 
                             ## although we can still explicitely request them 
    }
  } 
  ####
  #### REMLformula for iterativeFits is a false (by default full-model formula) REML formula for both iterative fits; that makes them *somewhat* more suitable for LR tests
  if (iterativeFits) {
    ## check fullm.list$REMLformula, which will be copied into nullm in all cases of fixed LRTs
    if (dotlist$HLmethod %in% c("ML","PQL/L")) {
      fullm.list$REMLformula <- NULL
      nullm.list$REMLformula <- NULL
    } else { ## an REML variant
      if (is.null(REMLformula)) { ## default
        fullm.list$REMLformula <- NULL ## so far this is the HLfit default
        nullm.list$REMLformula <- fullm.list$formula ## so far this is the HLfit default
      } else { ## allows alternative choice, but still the same *both* the full and the null fit
        fullm.list$REMLformula <- REMLformula
        nullm.list$REMLformula <- REMLformula
      }
      namesinit <- names(fullm.list$init.corrHLfit)
      namesinit <- setdiff(namesinit,c("rho","nu","Nugget","ARphi"))
      len <- length(namesinit)
      if ( len>0 && optimFits) {
        if (len > 1) namesinit <- paste(c(paste(namesinit[-len],collapse=", "),namesinit[len]),collapse=" and ")
        message("Argument 'init.corrHLfit' is used in such a way that")
        message(paste("  ",namesinit," will be estimated by maximization of p_v.",sep=""))
        message("  'REMLformula' will be inoperative if all dispersion")
        message("  and correlation parameters are estimated in this way.")
      }
    }
  }
  if (verbose["trace"]) {
    if (iterativeFits) {
      cat("corrMM.LRT will perform iterative fits of dispersion parameters with REML formula ")
      print(fullm.list$REMLformula)
    } else print("corrMM.LRT will not perform iterative fits of dispersion parameters.")  
    if (optimFits) {
      cat("corrMM.LRT will perform generic optimization fits of dispersion parameters with REML formula ")
      print(optimREMLformula)
    } else print("corrMM.LRT will not perform generic optimization fits of dispersion parameters.")
  }
  #########
  ## definitions for updating parameters
  ## local fn:
  `update.ranef.pars` <- function(from.fit,to.arglist,which.pars=c("rho","nu")) {
     ## default 'which.pars' prevents the dispersion params from being initialized by this fn
     ## using lambda leads to some poor results
     if ("lambda" %in% which.pars) {
        if ("lambda" %in% which.optim.fit ) {
          to.arglist$init.corrHLfit$lambda <- from.fit$lambda
        } else if ("lambda" %in% which.iterative.fit ) {
          to.arglist$init.HLfit$lambda <- from.fit$lambda
        } 
     } ## ELSE keeps any preexisting to.arglist$init.corrHLfit$lambda
     if ("phi" %in% which.pars) {
        if ("phi" %in% which.optim.fit ) {
          to.arglist$init.corrHLfit$phi <- from.fit$phi
        } else if ("phi" %in% which.iterative.fit ) {
          to.arglist$init.HLfit$phi <- from.fit$phi
        }
     }
     if ("rho" %in% which.pars) {
       if ("rho" %in% which.optim.fit )  { ## always by ML whether user provided initial values or not
         to.arglist$init.corrHLfit$rho <- from.fit$corrPars$rho   
       }    
     }  
     if ("ARphi" %in% which.pars) {
       if ("ARphi" %in% which.optim.fit )  { ## always by ML whether user provided initial values or not
         to.arglist$init.corrHLfit$ARphi <- from.fit$corrPars$ARphi   
       }    
     }  
     if ("nu" %in% which.pars) {
        if ("nu" %in% which.optim.fit) { ## always by ML whether user provided initial values or not
          to.arglist$init.corrHLfit$nu <- from.fit$corrPars$nu   
        }
     }      
     if ("v_h" %in% which.pars) { ## it is a very bad idea to update with v_h by default. Cf Ln... (binary)
        to.arglist$init.HLfit$v_h <- from.fit$v_h ## added 04/12/12
     }
    return(to.arglist)
  }
  ## another lcoal fn. ! different default which.pars
  init.fixed.ranefpars <- function(from.fit,to.arglist,which.pars=c("rho","nu","lambda","phi")) {
     ## default 'which.pars' prevents the dispersion params from being initialized by this fn
     ## using lambda leads to some poor results
     if ("lambda" %in% which.pars) {
        if ("lambda" %in% which.optim.fit ) {
          to.arglist$ranFix$lambda <- from.fit$lambda
        } else if ("lambda" %in% which.iterative.fit ) {
          to.arglist$ranFix$lambda <- from.fit$lambda
        } 
        to.arglist$init.corrHLfit$lambda <- NULL
     } ## ELSE keeps any preexisting to.arglist$init.corrHLfit$lambda
     if ("phi" %in% which.pars) {
        if ("phi" %in% which.optim.fit ) {
          to.arglist$ranFix$phi <- from.fit$phi
        } else if ("phi" %in% which.iterative.fit ) {
          to.arglist$ranFix$phi <- from.fit$phi
        }
        to.arglist$init.corrHLfit$phi <- NULL
     }
     if ("rho" %in% which.pars) {
        if ("rho" %in% which.optim.fit )  { ## always by ML whether user provided initial values or not
          to.arglist$ranFix$rho <- from.fit$corrPars$rho   
          to.arglist$init.corrHLfit$rho <- NULL
        }    
     }  
     if ("nu" %in% which.pars) {
       if ("nu" %in% which.optim.fit) { ## always by ML whether user provided initial values or not
         to.arglist$ranFix$nu <- from.fit$corrPars$nu   
         to.arglist$init.corrHLfit$nu <- NULL
       }
     }      
     if ("Nugget" %in% which.pars) {
       if ("Nugget" %in% which.optim.fit) { ## always by ML whether user provided initial values or not
         to.arglist$ranFix$Nugget <- from.fit$corrPars$Nugget   
         to.arglist$init.corrHLfit$Nugget <- NULL
       }
     }      
     return(to.arglist)
  }  
  ####
  #### ALWAYS p_v for testing fixed effects, p_bv for random effects 
  if ( ! is.null(null.predictor)) { ## ie if test effet fixe
    testFix <- T
    if (is.null(test.obj)) test.obj <- "p_v"
    if (! inherits(null.predictor,"predictor")) null.predictor <- Predictor(null.predictor)
    nullm.list$formula <- null.predictor
  } else if ( length(null.disp)>0 ) { ## test disp/corr param
    testFix <- F
    if (is.null(test.obj)) test.obj <- "p_bv"
    #
    mess <- pastefrom("changes to objective fn required to renullfit or refull fit computation.")
    stop(mess)
    #
    namescheck <- names(fullm.list$lower)
    namescheck <- namescheck[ ! namescheck %in% names(null.disp)]  
    nullm.list$lower <- nullm.list$lower[namescheck]
    nullm.list$upper <- nullm.list$upper[namescheck]
    nullm.list$ranFix <- c(nullm.list$ranFix,null.disp) ## adds fixed values to any preexisting one
  } else testFix <- NA
  ####
  #### run iterative fit
  ## FIRST NULL FIT; iterative then optim
  if (iterativeFits) { ## user requested not iterative (default), but in this case the algo still tries something for iterative lambda 
    nullREML.list <- nullm.list
    if (nullREML.list$HLmethod=="ML") { nullREML.list$HLmethod <- "REML" ## but do not over write more explicit HL(.,.) statements 
    } else if (substr(nullREML.list$HLmethod,0,2) == "ML") { nullREML.list$HLmethod <- paste("RE",substring(nullREML.list$HLmethod,3),sep="") ##  
    } else if (nullREML.list$HLmethod %in% c("PQL/L")) nullREML.list$HLmethod <- "REPQL" ## but do not over write more explicit HL(.,.) statements 
    nullREML.list$init.corrHLfit$lambda <- NULL 
    nullREML.list$lower$lambda <- NULL ## 2015/09/14 
    nullREML.list$upper$lambda <- NULL ## 2015/09/14 
    nullREML.list$init.corrHLfit$phi <- NULL ## 21/02/2013
    nullREML.list$control.HLfit$conv.threshold <- 1e-04 
    if (trace) nullREML.list$trace <- list(file="trace.lamREMLnullfit.txt",append=F)
    lamREMLnullfit <- do.call(method,nullREML.list) ## must return something that simulate() can manipulate
    trace.info <- data.frame(iter=0,step="lamREMLnullfit",obj=lamREMLnullfit$APHLs[[test.obj]])
    ## for ML 
    nullm.list <- update.ranef.pars(from.fit=lamREMLnullfit,to.arglist=nullm.list) 
    nullm.list <- update.ranef.pars(from.fit=lamREMLnullfit,to.arglist=nullm.list,which.pars=c("lambda")) 
  }  else lamREMLnullfit <- NULL 
  ####
  #### if not RE<...>, run optim fit
  if (optimFits) {
    if ("phi" %in% which.iterative) nullm.list$init.corrHLfit$phi <- NULL
    if ("lambda" %in% which.iterative) nullm.list$init.corrHLfit$lambda <- NULL
    nullm.list$REMLformula <- optimREMLformula ## ie ML fit... always the case anyway but needed if iterative...
    if (trace) nullm.list$trace <- list(file="trace.nullfit.txt",append=F)
    nullfit <- do.call(method,nullm.list) ## must return something that simulate() can manipulate
    if ( ! is.null(lamREMLnullfit)) { ## if an iterative fit was performed, we check whether it provided a better fit
        ## since lamREMLnullfit use full model formula and nullfit uses ML, one should not compare their p_v
        ## instead we fix the ranPars [since HLCor, not corrHLfit, call] and refit by ML
        MLforREMLranPars <- nullREML.list
        MLforREMLranPars$ranFix <- NULL
        MLforREMLranPars$REMLformula <- optimREMLformula
        MLforREMLranPars$ranPars<- with(lamREMLnullfit,c(corrPars,list(phi=phi,lambda=lambda)))
        MLforREMLranPars$coordinates <- dotlist$coordinates
        gnrf <- do.call("HLCor",MLforREMLranPars)$APHLs[[test.obj]]
        if (gnrf > nullfit$APHLs[[test.obj]]) { ## annoying as it means ML has failed to maximize p_v, but can occur
          trace.info <- rbind(trace.info,data.frame(iter=0,step="(!) nullfit (-)",obj=nullfit$APHLs[[test.obj]]))
          nullfit <-lamREMLnullfit ## step back
          ## this suggests a divergence of the inner loop; cf notes for 29/12/12 for case Hn
          ## Using which.iterative may be a fix
        } else trace.info <- rbind(trace.info,data.frame(iter=0,step="nullfit (+)",obj=nullfit$APHLs[[test.obj]]))
    } else trace.info <- rbind(trace.info,data.frame(iter=0,step="nullfit",obj=nullfit$APHLs[[test.obj]]))
  } else {
    nullfit <- lamREMLnullfit ### 
  }
## THEN FIRST FULL FIT, iterative then optim
  ## iterative
  fullm.list <- update.ranef.pars(from.fit=nullfit,to.arglist=fullm.list) 
  if (iterativeFits ) { ## in which case still tries somethng for iterative lambda 
    REML.list <- fullm.list
    if (REML.list$HLmethod=="ML") REML.list$HLmethod <- "REML" ## but do not over write more explicit HL(.,.) statements 
    if (REML.list$HLmethod %in% c("PQL/L")) REML.list$HLmethod <- "REPQL" ## but do not over write more explicit HL(.,.) statements 
    REML.list$init.corrHLfit$lambda <- NULL 
    REML.list$lower$lambda <- NULL ## 2015/09/14 
    REML.list$upper$lambda <- NULL ## 2015/09/14 
    REML.list$init.corrHLfit$phi <- NULL ## 21/02/2013
    REML.list$control.HLfit$conv.threshold <- 1e-04 
    if (trace) REML.list$trace <- list(file="trace.lamREMLfullfit.txt",append=F)
    lamREMLfullfit <- do.call(method,REML.list) ## this performs an REML fit for lambda, provide the initial value for the ML fit
    trace.info <- rbind(trace.info,data.frame(iter=0,step="lamREMLfullfit",obj=lamREMLfullfit$APHLs[[test.obj]]))
    ## for ML
    fullm.list <- update.ranef.pars(from.fit=lamREMLfullfit,to.arglist=fullm.list)  
    fullm.list <- update.ranef.pars(from.fit=lamREMLfullfit,to.arglist=fullm.list,which.pars=c("lambda")) 
  } else lamREMLfullfit <- NULL 
  if (optimFits) {
    if ("phi" %in% which.iterative) fullm.list$init.corrHLfit$phi <- NULL
    if ("lambda" %in% which.iterative) fullm.list$init.corrHLfit$lambda <- NULL
    fullm.list$REMLformula <- optimREMLformula ## ie ML fit... always the case anyway but needed if iterative...
    if (trace) fullm.list$trace <- list(file="trace.fullfit.txt",append=F)
    fullfit <- do.call(method,fullm.list) ## must return something that simulate() can manipulate
    if ( ! is.null(lamREMLfullfit)) { ## if an iterative fit was performed, we check whether it provided a better fit
        MLforREMLranPars <- REML.list
        MLforREMLranPars$ranFix <- NULL
        MLforREMLranPars$REMLformula <- optimREMLformula
        MLforREMLranPars$ranPars<- with(lamREMLfullfit,c(corrPars,list(phi=phi,lambda=lambda)))
        MLforREMLranPars$coordinates <- dotlist$coordinates
        gnrf <- do.call("HLCor",MLforREMLranPars)$APHLs[[test.obj]]
        if (gnrf > fullfit$APHLs[[test.obj]]) { ## annoying as it means ML has failed to maximize p_v, but can occur
          trace.info <- rbind(trace.info,data.frame(iter=0,step="(!) fullfit (-)",obj=fullfit$APHLs[[test.obj]]))
          fullfit <-lamREMLfullfit
        } else trace.info <- rbind(trace.info,data.frame(iter=0,step="fullfit (+)",obj=fullfit$APHLs[[test.obj]]))
    } else trace.info <- rbind(trace.info,data.frame(iter=0,step="fullfit",obj=fullfit$APHLs[[test.obj]]))
  } else {
    fullfit <- lamREMLfullfit
  }
  ## ITERATIONS
  nullnamesX <- colnames(nullfit$`X.pv`)
  fullnamesX <- colnames(fullfit$`X.pv`)
  nullVarsPos <- which( fullnamesX %in% nullnamesX) 
  profileVars <- which( ! fullnamesX %in% nullnamesX) ## les parametres testes
  profileX <- fullfit$`X.pv`[,profileVars,drop=F]
  locit <- 1
  if (! (optimFits && locit <= maxit) ) { ## will not enter the loop => LRTori ici
    LRTori <- 2*(fullfit$APHLs[[test.obj]]-nullfit$APHLs[[test.obj]])
    conv <- 0 ## will skip the loop
  } else  conv <- 1e08 ## LRTori will be computed at first iteration of the loop
  newnull <- T; newfull <- T; 
  newforprof <- T
  LRTprof <- NA
  ## must *enter* the loop *iff* (optimFits)
  while ( optimFits && conv>0.001 && (newnull || newfull) && locit <= maxit ) { ## optimFits required since the later fits use optimREMLformula; but the premise could be reconsidered
    ####### this loop aims to ensure that maximal nullfit and fullfit are found
    ## currently implemented only in cases where some random effect parameters are fit
    ## (test below; no alternative yet)
    #######
    conv <-0 
    prof.needed <- F ## must be reset at the beginning of each loop
    newnull <- F
    newfull <- F
    if (length(c(which.iterative.fit,which.optim.fit))>0) { ## some random effect parameters are fit 
      ####### while fullfit or nullfit has been changed in the last iter
      # if (fullfit seems OK) then use it to refine nullfit: {
      #   renullfit with update.ranef.pars
      #   renullfit with further update lambda, and iter.mean.dispFix <- 10
      #   retain the best (so that update lambda, and iter.mean.dispFix <- 10 can contaminate later renullfit and the profile) 
      # }
      # if (fullfit lower than refined nullfit) then use nullfit to refine fullfit: {
      #   refullfit with update.ranef.pars and init.HLfit$fixef[nullVarsPos] <- as.numeric(nullfit$fixef)
      #   renullfit with further update lambda, and iter.mean.dispFix <- 10
      #   retain the best (so that update lambda, and iter.mean.dispFix <- 10 can contaminate later renullfit) 
      # }
      #######      
if (restarts) { 
      if ( fullfit$APHLs[[test.obj]] > nullfit$APHLs[[test.obj]]) { 
        ## ideally always the case. Then we look whether we can improve nullfit using fullfit as starting point
        ## (else, we will improve fullfit first)
        ## REFIT NULL FIT
        renullm.list <- nullm.list
        if (trace) renullm.list$trace <- list(file="trace.renullfit.txt",append=F)
        renullm.list <- update.ranef.pars(from.fit=fullfit,to.arglist=renullm.list)
        if ("phi" %in% which.iterative) renullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) renullm.list$init.corrHLfit$lambda <- NULL
        renullfit <- do.call(method,renullm.list) 
        if (renullfit$APHLs[[test.obj]] > nullfit$APHLs[[test.obj]]) {
          conv <- conv + renullfit$APHLs[[test.obj]] - nullfit$APHLs[[test.obj]]
          nullfit <-renullfit
          nullm.list <- renullm.list ## updated for immediate use and use by profilize
          newnull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="renullfit.1 (+)",obj=nullfit$APHLs[[test.obj]]))
        } 
        if (trace) renullm.list$trace <- list(file="trace.renullfit.txt",append=T)
        renullm.list <- update.ranef.pars(from.fit=fullfit,to.arglist=renullm.list)
        renullm.list <- update.ranef.pars(from.fit=fullfit,to.arglist=renullm.list,which.pars=c("lambda"))
        ############################ renullm.list$control.HLfit$iter.mean.dispFix <- 10
        if ("phi" %in% which.iterative) renullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) renullm.list$init.corrHLfit$lambda <- NULL
        renullfit <- do.call(method,renullm.list) 
        if (renullfit$APHLs[[test.obj]] > nullfit$APHLs[[test.obj]]) {
          conv <- conv + renullfit$APHLs[[test.obj]] - nullfit$APHLs[[test.obj]]
          nullfit <-renullfit
          nullm.list <- renullm.list ## updated for use by profilize
          newnull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="renullfit.2 (+)",obj=nullfit$APHLs[[test.obj]]))
        } 
      } 
      ## two heuristic tries for large improvements, else then a sure small improvement
      if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]) { ## ideally never the case
        ## REFIT FULL FIT
        ## then we try other starting values
        ## the starting lambda value seems to be quite important
        refullm.list <- fullm.list
        if (trace) refullm.list$trace <- list(file="trace.refullfit.txt",append=F)
        refullm.list <- update.ranef.pars(from.fit=nullfit,to.arglist=refullm.list)
        refullm.list$init.HLfit$fixef <- rep(0,length(fullnamesX)) ## fixef must be modified here
#        refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(nullfit$fixef) 
        refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(mean(nullfit$eta)) 
        if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
        refullfit <- do.call(method,refullm.list)
        if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## 
          conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]] ## only D(FULL)
          fullfit <- refullfit
          fullm.list <- refullm.list ## updated for immediate use 
          newfull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="refullfit.1 (+)",obj=fullfit$APHLs[[test.obj]]))
        }
      }
      if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]) { ## if previous heuristic attempt not good enough
        ####### try something else
        if (trace) refullm.list$trace <- list(file="trace.refullfit.txt",append=T)
        refullm.list <- update.ranef.pars(from.fit=nullfit,to.arglist=refullm.list,which.pars=c("lambda"))
        ## ie, essentially, refullm.list$init.corrHLfit$lambda <- nullfit$lambda, which works *sometimes*          
        ################## refullm.list$control.HLfit$iter.mean.dispFix <- 10
        if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
        refullfit <- do.call(method,refullm.list)
        if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## 
          conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]] ## only D(FULL)
          fullfit <- refullfit
          fullm.list <- refullm.list ## maybe not useful
          newfull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="refullfit.2 (+)",obj=fullfit$APHLs[[test.obj]]))
        } 
      }
      if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]] ) { ## if previous heuristic attemptS not good enough
        ## Formally a call to corrHLfit, but all optim parameters are fixed. Hence this is only an HLCor call, and it generates no valid output in the trace file
        fullm.con.list <- fullm.list
        fullm.con.list <- init.fixed.ranefpars(from.fit=nullfit,to.arglist=fullm.con.list)
        fullm.con.list <- update.ranef.pars(from.fit=nullfit,to.arglist=fullm.con.list,which.pars=c("v_h"))
        fullm.con.list$init.HLfit$fixef <- rep(0,length(fullnamesX))
#        fullm.con.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(nullfit$fixef) 
        fullm.con.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(mean(nullfit$eta)) 
        if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
        if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
        if (trace) fullm.con.list$trace <- list(file="trace.refullfit.txt",append=T)
        refullfit <- do.call(method,fullm.con.list)
        if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## but not necess better than nullfit
          conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]] ## only D(FULL)
          fullfit <- refullfit
          fullm.list <- update.ranef.pars(from.fit=refullfit,to.arglist=fullm.list)
          fullm.list <- update.ranef.pars(from.fit=refullfit,to.arglist=fullm.list,which.pars=c("lambda"))
          fullm.list$init.HLfit$fixef <- as.numeric(refullfit$fixef) 
          newfull <- T
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="refullfit.3 (+)",obj=fullfit$APHLs[[test.obj]]))
        } 
        if ( fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]-0.001) { ## if still BAD
          ## serious problem with iterative algo in the simplest case
          ## profiling is our last chance
          prof.needed <- T ## => fullfit$APHLs[[test.obj]] < nullfit$APHLs[[test.obj]]-0.001 
          trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.needed",obj=fullfit$APHLs[[test.obj]]))
        }
      }
}     
      if (locit==1L) LRTori <- 2*(fullfit$APHLs[[test.obj]]-nullfit$APHLs[[test.obj]]) ## first iteration only       
      if (newfull || newnull) {
        newforprof <- T ## T at initialization and after any 'unused' new full/null
      } else if (locit>1) trace.info <- rbind(trace.info,data.frame(iter=locit,step="no.new(f/n)ull",obj=fullfit$APHLs[[test.obj]]))
      if ( ! newforprof) {
        if (prof.needed) {
          print(trace.info)
          mess <- pastefrom("prof.needed && ! (newnull || newfull).")
          message(mess)
          ## conv must be null (else one the 'new...' is true) hence we will exit the loop...
        }
      } else if ( prof.needed || profiles>0 ) { ## given newforprof
        if ( ! prof.needed ) profiles <- profiles - 1 ## 'pays' one profile
        # then we need to maximize p_v(tau(beta)) over different beta values 
        # and therefore we need to be able to fit for fixed beta values 
        ## il faut (1) identifier les fixef sur lesquel le profile est fait
        prof.list <- nullm.list ## uses last null fit here
        if (trace) {
          try(unlink("trace.profile.txt")) ## the file is written in by HLCor()                   
          prof.list$trace <- list(file="trace.profile.txt",append=T)
        }
        profilize <- function(offsetbeta) {
          addOffset <- profileX %*% offsetbeta ## semble marcher pour scalaire comme pour vecteur
          offset <- attr(prof.list$formula,"offsetObj")$total
          if ( is.null(offset)) stop("problem with offset in profiling code")
          attr(prof.list$formula,"offsetObj") <- list(total=offset + addOffset,nonZeroInfo=TRUE)
          prof.fit <- do.call(method,prof.list) ## 
          return(prof.fit$APHLs[[test.obj]])
        }
        beta_se <- sqrt(diag(fullfit$beta_cov))
        if (length(profileVars)>1) {
          ## optim syntax
          ## uses last full fit here:
          # lowup <- fullfit$fixef[profileVars] %*% t(c(-1,5)/2) ## should be of dim [,2]
          # lowup <- apply(lowup,1,sort) ##  ncol=2
          lowup <- fullfit$fixef[profileVars] + beta_se[profileVars] %*% t(c(-2,2)) ## ncol=2
          lowp <- lowup[,1]; upp <- lowup[,2]
          ## if we are here because nullfit > fullfit, nullfit should be in the interval  
          lowp[lowp>0] <- - (beta_se[profileVars])[lowp>0]
          upp[upp<0] <- (beta_se[profileVars])[upp<0]
          stop("optim call missing 'here' in corrMM.LRT")
          profmax <- optim(par=fullfit$fixef[profileVars],
                              fn=profilize,lower=lowp,upper=upp,
                              control=list(fnscale=-1),method="L-BFGS-B")
          mess <- pastefrom("code needed here for optim output.")
          stop(mess)
        } else {
          ## one dim optimize syntax
## ad hoc code of test versions <= 1.0
#          if (fullfit$fixef[profileVars] != 0) { ## if fullfit sufficiently different from null fit
#            interval <- sort(fullfit$fixef[profileVars]*c(-1,5)/2)
#          } else interval <- 0.001*c(-1,1) ## very quick ad hoc patch
          interval <- fullfit$fixef[profileVars]+c(-2,2)*beta_se[profileVars]   
          ## if we are here because nullfit > fullfit, nullfit should be in the interval  
          if (interval[1]>0) interval[1] <- - beta_se[profileVars]
          if (interval[2]<0) interval[2] <- beta_se[profileVars]
          profmax <- optimize(profilize,interval=interval,maximum=T)  ## FR->FR optim avec un valuer init inteeligente ?
          addOffset <- profileX %*% profmax$maximum
          offset <- attr(prof.list$formula,"offsetObj")$offset
          if ( is.null(offset)) stop("problem with offset in profiling code") ## ? FR->FR code repetitif, creer update.offset
          attr(prof.list$formula,"offsetObj") <- list(total=offset + addOffset,nonZeroInfo=TRUE)
          proffit <- do.call(method,prof.list) ## to recover the ranef parameters etc
          if ( proffit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## if profile max should improve fullfit
            trace.info <- rbind(trace.info,data.frame(iter=locit,step="proffit (+)",obj=proffit$APHLs[[test.obj]]))
            conv <- conv + proffit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]]
            refullm.list <- update.ranef.pars(from.fit=proffit,to.arglist=fullm.list)
            refullm.list <- update.ranef.pars(from.fit=proffit,to.arglist=refullm.list,which.pars=c("lambda","phi"))
            ## apparently no a good idea to provide init.HLfit$fixef without $v_h
#            refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(proffit$fixef) 
            refullm.list$init.HLfit$fixef[nullVarsPos] <- as.numeric(mean(proffit$eta)) 
            refullm.list$init.HLfit$fixef[profileVars] <- as.numeric(profmax$maximum) 
            refullm.list$init.HLfit$v_h <- proffit$v_h 
            if ("phi" %in% which.iterative) refullm.list$init.corrHLfit$phi <- NULL
            if ("lambda" %in% which.iterative) refullm.list$init.corrHLfit$lambda <- NULL
            refullfit <- do.call(method,refullm.list) 
            if ( refullfit$APHLs[[test.obj]] > fullfit$APHLs[[test.obj]] ) { ## if confirmed improvement of fullfit 
              conv <- conv + refullfit$APHLs[[test.obj]] - fullfit$APHLs[[test.obj]]
              fullfit <- refullfit
              fullm.list <- refullm.list ## maybe not useful
              newfull <- T
              trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.refullfit (+)",obj=fullfit$APHLs[[test.obj]]))
            } else if ( refullfit$APHLs[[test.obj]] < proffit$APHLs[[test.obj]]-0.001 ) { ## if failure to confirm improvement of fullfit
              trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.refullfit (-)",obj=refullfit$APHLs[[test.obj]]))
              mess <- pastefrom("refullfit not as good as profile fit.")
              message(mess)
            } else { ## refull < full < proffit but all within 0.001 logL units 
              trace.info <- rbind(trace.info,data.frame(iter=locit,step="prof.refullfit (=)",obj=refullfit$APHLs[[test.obj]]))
            } 
          } else if ( proffit$APHLs[[test.obj]] < fullfit$APHLs[[test.obj]]-0.001 ) { ## somewhat problematic
            trace.info <- rbind(trace.info,data.frame(iter=locit,step="proffit (-)",obj=proffit$APHLs[[test.obj]]))
          } else {
            trace.info <- rbind(trace.info,data.frame(iter=locit,step="proffit (=)",obj=proffit$APHLs[[test.obj]]))
          }
        }
        newforprof <- F ## last newfull/newnull have been 'used'
        LRTprof <- 2*(fullfit$APHLs[[test.obj]]-nullfit$APHLs[[test.obj]])
      }
    } ## dispersion/correlation params were fit
    locit <- locit + 1
  } ## end loop
  ## BOOTSTRAP
  if ( ! is.na(testFix)) {
    if (boot.repl>0) {
      ## the following prevents a bug when we change "this" function name (but works only for certain calls...)
      thisFnName <- as.character(sys.call()[[1]]) 
      bootlist <- dotlist ## copies ranFix
      bootlist$control <- bootcontrol ## (a list)
      bootlist <- c(bootlist,list(null.predictor=null.predictor,null.disp=null.disp,REMLformula=REMLformula,method=method)) ## unchanged user REMLformula forwarded
      bootlist$verbose <- c(summary=FALSE)
      bootlist$trace <- FALSE 
      bootlist$boot.repl <- 0 ## avoids recursive call of bootstrap
      ## all.optim.vars <- c(which.optim.fit,which.iterative.fit) ## looks suspect, but is correct because if there is no dispersion estimate in  
      ## init.corrHLfit from the corrMM call, then optimFits has been set to F, 
      ## and then we perform only iterativeFits which takes out the dispersion estimates from init.corrHLfit.
      ## A more transparent code is:
      if (optimFits) {
        all.optim.vars <- c(which.optim.fit,which.iterative.fit)
      } else all.optim.vars <- which.optim.fit
      for (st in all.optim.vars) {
        if ( st %in% names(nullfit$corrPars)) {
          bootlist$init.corrHLfit[st] <- nullfit$corrPars[st]
        } else if ( st %in% names(nullfit)) {
          bootlist$init.corrHLfit[st] <- nullfit[st] ## handled in dotlist by corrHLfit, cf notes 090113
        } ## it's also possible that st is nowhere (10/2013: currently for Nugget)
      }
      bootreps<-matrix(,nrow=boot.repl,ncol=2) 
      colnames(bootreps) <- paste(c("full.","null."),test.obj,sep="")
      cat("bootstrap replicates: ")
      simbData <- nullfit$data
      if (tolower(nullfit$family$family)=="binomial") {
        form <- attr(nullfit$predictor,"oriFormula") ## this must exists...  
        if (is.null(form)) {
          mess <- pastefrom("a 'predictor' object must have an 'oriFormula' member.",prefix="(!) From ")
          stop(mess)
        }
      }
      ## the data contain any original variable not further used; e.g original random effect values in the simulation tests  
      for (ii in 1:boot.repl) {
        locitError <- 0
        repeat { ## for each ii!
          newy <- simulate(nullfit) ## cannot simulate all samples in one block since some may not be analyzable  
          if (tolower(nullfit$family$family)=="binomial") {
            ## We have different possible cbind(exprL,exprR) arguments, but in all case the predictor is that of exprL and 
            #  exprR is $weights- exprL 
            
            ## c'est bouseux: soit j'ai (pos, neg) et le remplacement est possible
            ##    soit j'ai (pos,ntot -pos) et le 2e remplacment n'est pas poss (et pas necess)
            ##    aussi (ntot - pos, pos) ...
            ## would be simple if always ntot-pos, but how to control this ? 
            ## simbData[[as.character(form[[2]][[2]])]] <- newy
            ## simbData[[as.character(form[[2]][[3]])]] <- nullfit$weights - newy    
            exprL <- as.character(form[[2]][[2]]) 
            exprR <- as.character(form[[2]][[3]]) 
            if (length(exprL)==1L) simbData[[exprL]] <- newy 
            if (length(exprR)==1L) simbData[[exprR]] <- nullfit$weights - newy                    
            ## if (length(exprR)! =1) exprRdoes not correspond to a column in the data.frame so there is no column to replace                     
          } else {simbData[[as.character(nullfit$predictor[[2]])]] <- newy}
          bootlist$data <- simbData
          bootrepl <- try(do.call(thisFnName,bootlist))  ################# CALL ################# 
          if ( ! inherits(bootrepl,"try-error")) { ## eg separation in binomial models... alternatively, test it here (require full and null X.pv... )
            bootreps[ii,] <- c(bootrepl$fullfit$APHLs[[test.obj]],bootrepl$nullfit$APHLs[[test.obj]])
            break ## replicate performed, breaks the repeat
          } else { ## there was one error
            locitError <- locitError + 1
            if (locitError>10) { ## to avoid an infinite loop
              stop("Analysis of bootstrap samples fails repeatedly. Maybe no statistical information in them ?")
            } ## otherwise repeat!
          }
        } 
        cat(ii);cat(" ")
        if ((ii %% 50)==0L) cat("\n")
      } ## end main bootstrap loop
      cat("\n") ##  
    } ## end bootstrap 
  } else { ## nothing operativ yet
    bootreps<-matrix(,nrow=boot.repl,ncol=length(unlist(fullfit$APHLs))) 
    colnames(bootreps) <- names(unlist(fullfit$APHLs))
    ## more needed here ?
  }
  ## prepare output (messy for back compatibility)
  if ( ! is.na(testFix)) {
    if (testFix) {df <- length(fullfit$fixef)-length(nullfit$fixef)} else {df <- length(null.disp)}
    resu <- list(fullfit=fullfit,nullfit=nullfit)
    LRTinfo <- list(df=df,LRTprof = LRTprof,LRTori = LRTori)
    if (boot.repl>0) {
      meanbootLRT <- 2*mean(bootreps[,1]-bootreps[,2]) 
      LRTinfo$meanbootLRT <- meanbootLRT
      LRTinfo$bootreps <- bootreps
    }
    basicLRT <- data.frame(LR2=LRTori,df=df,pvalue=1-pchisq(LRTori,df=df))
  } else {
    resu <- list(fullfit=fullfit)
    LRTinfo <- list()
    basicLRT <- list()
  }
  LRTinfo$trace.info <- trace.info 
##  resu$LRTinfo <- LRTinfo ## pas compatible avec hglmjob.R...
  resu <- c(resu,LRTinfo) ## loses the sublist structure, which wouldnot be compatible with hglmjob.R...  
  resu <- c(resu,list(basicLRT = basicLRT)) ## as in anova and LRT and their summaries  
  class(resu) <- c("fixedLRT",class(resu)) 
  return(resu)
}

summary.fixedLRT <- function(object,verbose=TRUE,...) {
  LRT <- object$basicLRT
  print(object$basicLRT)
  #   if (verbose) {
  #     cat(" ========      'full' model:     ========\n")    
  #     summary(object$fullfit,...) 
  #     cat(" ========      'null' model:     ========\n")    
  #     summary(object$nullfit,...) 
  #     cat(" ======== Likelihood ratio test: ========\n")    
  #   }
  #   outst <- paste(" LR statistic (",object$df," df): ",signif(object$LRTori,3),sep="")    
  #   cat(outst)
  bootInfo <- object$bootInfo
  if (!is.null(bootInfo)) {
    cat(" ======== Bootstrap: ========\n")    
    outst <- paste("Raw simulated p-value: ",signif(object$rawBootLRT$pvalue,3),sep="")    
    cat(outst)
    X2 <- object$BartBootLRT$LR2
    outst <- paste("\nBartlett-corrected LR statistic (",object$BartBootLRT$df," df): ",signif(X2,3),sep="")    
    cat(outst)
  } 
}

print.fixedLRT <-function(x,...) {
  summary(x,...)
  invisible(x)
}

    
