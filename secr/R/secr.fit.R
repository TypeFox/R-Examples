################################################################################
## package 'secr'
## secr.fit.R
## moved from methods.R 2011-01-30
## 2011-10-20 generalized designD
## 2011-10-20 renamed 'grps' 'grouplevels'
## 2011-12-16 streamlined preparation of models
## 2012-01-22 purged phi/turnover
## 2012-01-31 experimental addition of parameter cut
## 2012-04-06 'fixed' bug fixed (see functions.r)
## 2012-07-24 unmash component of details
## 2013-04-11 hcov
## 2013-04-19 lambda0
## 2013-04-20 default mask type changed to trapbuffer
## 2013-07-02 esa parameterisation details$param = 2
## 2013-07-19 a0 parameterisation details$param = 3
## 2013-10-14 tweak a0 biasD
## 2013-10-28 revise pmix models
## 2013-11-09 rearrange model check code to put in warning
## 2014-03-12 sigmak parameterisation details$param = 4
## 2014-03-12 bufferbiascheck now in suggest.buffer
## 2014-03-24 allow combination esa, sigmak (param 6)
## 2014-06-02 secr.design.MS argument ignoreusage
## 2014-08-22 smoothers... badsmooth check
## 2014-08-25 model validation in separate function valid.model (utility.r)
## 2014-08-27 default model = NULL (model = list(D~1, g0~1, sigma~1)) REVERSED 2014-09-23
## 2014-09-06 warning if use Euclidean distances with linearmask
## 2014-09-23 directly replace g0 by lambda0 for detectfn 14:18
## 2014-10-13 designNE for non-Euclidean parameter
## 2014-10-14 allow start to be incomplete list of real parameter values
## 2014-10-16 rewrite of pnames and model section for clarity
## 2014-12-04 tweak to avoid RPSV problem with zero captures
## 2015-01-06 catch 'too few detections' before autoini
## 2015-04-01 details$autoini
## 2015-04-12 default telemetrytype 'concurrent' was overwritten 'none'
## 2015-05-24 default details$minprob changed to 1e-200
## 2015-10-09 reinstated pID and mark-sight
## 2015-10-11 adjusted default start values polygon detectors
###############################################################################

  secr.fit <- function (capthist,  model = list(D~1, g0~1, sigma~1), mask = NULL,
    buffer = NULL, CL = FALSE, detectfn = NULL, binomN = NULL, start = NULL,
    link = list(), fixed = list(), timecov = NULL, sessioncov = NULL, hcov = NULL,
    groups = NULL, dframe = NULL, details = list(), method = 'Newton-Raphson',
    verify = TRUE, biasLimit = 0.01, trace = NULL, ncores = 1, ...)

{
# Fit spatially explicit capture recapture model
#
# Arguments:
#
#  capthist   -  capture history object (includes traps object as an attribute)
#  model      -  formulae for real parameters in terms of effects and covariates
#  mask       -  habitat mask object
#  buffer     -  default buffer width should mask not be provided
#  CL         -  logical switch : conditional likelihood (T) or full likelihood (F)
#  detectfn   -  code for detection function 0 = halfnormal, 1 = hazard, 2 = exponential etc.
#  start      -  start values for maximization (numeric vector link scale);
#                if NULL then 'autoini' function is used
#  link       -  list of parameter-specific link function names 'log', 'logit', 'identity',
#                'sin', 'neglog'
#  fixed      -  list of fixed values for named parameters
#  timecov    -  data for time covariates if these are used in 'model'
#  sessioncov -  dataframe of session-level covariates
#  groups     -  vector of names to group fields in attr(capthist,'covariates') dataframe
#  dframe     -  optional data frame of design data for detection model (tricky & untested)
#  details    -  list with several additional settings, mostly of special interest
#  method     -  optimization method (indirectly chooses
#  verify     -  logical switch for pre-check of capthist and mask with verify()
#  trace      -  logical; if TRUE output each likelihood as it is calculated
#  ...        -  other arguments passed to nlm() or optim()


    #################################################
    ## Remember start time and call

    ptm  <- proc.time()
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")

    cl   <- match.call(expand.dots = TRUE)

    ## 2014-02-13
    if (is.character(capthist)) {
        capthist <- get(capthist, pos=-1)
    }
    if (is.character(mask)) {
        mask <- get(mask, pos=-1)
    }
    if (is.character(dframe)) {
        dframe <- get(dframe, pos=-1)
    }

    if (!inherits(capthist, 'capthist'))
        stop ("requires 'capthist' object")

    #################################################
    ## Default detection function

    if (is.null(detectfn)) {
        if (detector(traps(capthist)) %in% c('signal')) {
            detectfn <- 10
            warning ("detectfn not specified; using signal strength (10)")
        }
        else if (detector(traps(capthist)) %in% c('signalnoise')) {
            detectfn <- 12
            warning ("detectfn not specified; using signal-noise (12)")
        }
        else {
            detectfn <- 0
        }
    }

    else {
        if (detector(traps(capthist)) == 'presence')
            detectfn <- valid.detectfn(detectfn, 0:8)
        else
            detectfn <- valid.detectfn(detectfn)
    }

    #################################################
    ## Use input 'details' to override various defaults

    defaultdetails <- list(distribution = 'poisson',
                           hessian = 'auto',
                           trace = TRUE,
                           LLonly = FALSE,
                           centred = FALSE,
                           binomN = 1,
                           cutval = 0,
                           ## 2015-05-24 minprob = 1e-50,
                           minprob = 1e-200,
                           tx = 'identity',
                           param = 0,
                           unmash = FALSE,
                           telemetrytype = 'concurrent',
                           telemetrysigma = FALSE,
                           telemetrybvn = FALSE,
                           ignoreusage = FALSE,
                           debug = FALSE,
                           intwidth2 = 0.8,
                           normalize = FALSE,
                           usecov = NULL,
                           userdist = NULL,
                           autoini = 1,
                           knownmarks = TRUE,
                           nsim = 0,
                           chatonly = FALSE,
                           chat = NULL
                           )

    if (detector(traps(capthist)) %in% .localstuff$countdetectors)
        defaultdetails$binomN <- 0   ## Poisson
    if (!is.null(attr(capthist,'cutval')))
        defaultdetails$cutval <- attr(capthist,'cutval')
    else if (ms(capthist) & !is.null(attr(capthist[[1]],'cutval')))   ## 2012-09-04
        defaultdetails$cutval <- attr(capthist[[1]],'cutval')
    if (is.logical(details$hessian))
        details$hessian <- ifelse(details$hessian, 'auto', 'none')
## details$telemetrytype <- match.arg(details$telemetrytype, c('none', 'independent',
##                           'dependent', 'concurrent'))
## 2015-04-12, 2015-06-14
    details$telemetrytype <- match.arg(details$telemetrytype,
        c('concurrent', 'none', 'independent', 'dependent'))
    if (is.null(attr(capthist,'xylist')))
         details$telemetrytype <- 'none'
    if (details$telemetrytype == 'independent')
        details$telemetrysigma <- TRUE
    details <- replace (defaultdetails, names(details), details)
    if (details$telemetrytype == 'none' && details$telemetrysigma == TRUE)
        details$telemetrytype <- 'independent'
    if (!is.null(trace)) details$trace <- trace
    if (!is.null(binomN)) {
        if (detector(traps(capthist)) == 'count') {
            if (tolower(binomN) == 'usage')
                binomN <- 1   ## code for 'binomial size from usage' 2012-12-22
            if (tolower(binomN) == 'poissonhazard')
                binomN <- -1  ## interpret g() as true detection function 2013-01-07
        }
        details$binomN <- binomN   ## 2011 01 28
    }
    if (details$LLonly)  details$trace <- FALSE

    if (!(detector(traps(capthist)) %in% c('single','multi')) & (details$param == 1)) {
        warning ("Gardner & Royle parameterisation not appropriate, using param = 0")
        details$param <- 0
    }

    ## 2011-02-06 to be quite clear -
    if (detector(traps(capthist)) %in% c(.localstuff$exclusivedetectors,
                                         'proximity','signal','signalnoise'))
        details$binomN <- 1;

    #################################################
    ## MS - indicator TRUE if multi-session (logical)
    ## sessionlevels - names of sessions (character)

    MS <- ms(capthist)
    sessionlevels <- session(capthist)

    if (is.null(sessionlevels)) sessionlevels <- '1'
    anycount <- any(detector(traps(capthist)) %in% .localstuff$countdetectors)
    anypoly  <- any(detector(traps(capthist)) %in% c('polygon',  'polygonX'))
    anytrans <- any(detector(traps(capthist)) %in% c('transect', 'transectX'))
    alltelem <- all(detector(traps(capthist)) %in% c('telemetry'))
    if (alltelem) CL <- TRUE

    if (MS) {
       if (any (sapply(traps(capthist), detector) == 'single'))
        warning ("multi-catch likelihood used for single-catch traps")
    }
    else {
       if (detector(traps(capthist)) == 'single')
        warning ("multi-catch likelihood used for single-catch traps")
    }

    #################################################
    ## Optional data check added 2009 09 19

    if (verify) {
        memo ('Checking data', details$trace)
        test <- verify(capthist, report = 1)
        if (test$errors)
            stop ("'verify' found errors in 'capthist' argument")

        if (!is.null(mask)) {
            if (MS & ms(mask)) {
                ## list of masks
                test <- lapply(mask, verify, report = 1)
                notOK <- any(unlist(test))
            }
            else notOK <- verify(mask, report = 1)$errors
            if (notOK)
                stop ("'verify' found errors in 'mask' argument")
        }
    }

    #################################################
    ## Ensure valid mask
    ## assume traps(capthist) will extract a list of trap layouts
    ## if multi-session (MS == TRUE)

    usebuffer <- is.null(mask)    ## flag for later check
    if (usebuffer) {
        if (is.null(buffer)) {
            buffer <- 100
            if (!(detector(traps(capthist))=='presence') & !alltelem)
                warning ("using default buffer width 100 m")
        }
        if (MS) mask <- lapply (traps(capthist), make.mask, buffer = buffer, type = "trapbuffer")
        else    mask <- make.mask(traps(capthist), buffer = buffer, type = "trapbuffer")
    }
    else {
      if (MS & !ms(mask)) {
          if (inherits(mask, 'linearmask'))
              newclass <- c('list', 'linearmask', 'mask')
          else
              newclass <- c('list', 'mask')
          ## inefficiently replicate mask for each session!
          mask <- lapply(sessionlevels, function(x) mask)
          class (mask) <- newclass
          names(mask) <- sessionlevels
      }
    }

    nc <- ifelse (MS, sum(sapply(capthist, nrow)), nrow(capthist))
    if (nc < 1)
        warning (nc, " detection histories")

    if (is.null(details$userdist) & inherits(mask, 'linearmask'))
        warning ("using Euclidean distances with linear mask")

    #################################################
    ## mark-resight

    if (MS) {
        sighting <- sighting(traps(capthist[[1]]))
        Tu <- Tu(capthist[[1]])
        Tm <- Tm(capthist[[1]])
    }
    else {
        sighting <- sighting(traps(capthist))
        Tu <- Tu(capthist)
        Tm <- Tm(capthist)
    }
    if (('pID' %in% names(fixed)) & !is.null(Tm)){
        if ((fixed$pID == 1) & (sum(Tm)>0))
        warning ("mark-resight nonID sightings ignored when fixed$pID = 1")
    }
    if (sighting & CL & !is.null(Tu)) {
        warning ("mark-resight unmarked (but not nonID) sightings ignored when CL = TRUE")
    }

    #################################################
    ## optional centring of traps and mask 2010 04 27
    if (details$centred) {
        centre <- function (xy, dxy) {
            xy[,] <- sweep(xy, MARGIN = 2, FUN='-', STATS = dxy)
            xy
        }
        if (MS) {
            nsess <- length(traps(capthist))
            offsetxy <- lapply(traps(capthist), function(xy) apply(xy, 2, mean))
            for (i in 1:nsess) {
                temptraps <- centre(traps(capthist[[i]]), offsetxy[[i]])
                traps(capthist[[i]]) <- temptraps
                mask[[i]] <- centre(mask[[i]], offsetxy[[i]])
                attr(mask[[i]], 'meanSD')[1,1:2] <- attr(mask[[i]], 'meanSD')[1,1:2] -
                    offsetxy[[i]]
                attr(mask[[i]], 'boundingbox') <- centre(attr(mask[[i]], 'boundingbox'),
                    offsetxy[[i]])
            }
        }
        else {
            offsetxy <- apply(traps(capthist), 2, mean)
            traps(capthist) <- shift(traps(capthist), -offsetxy)
            mask <- shift.traps(mask, -offsetxy)
            attr(mask, 'meanSD')[1,1:2] <- attr(mask, 'meanSD')[1,1:2] - offsetxy
            attr(mask, 'boundingbox') <- centre(attr(mask, 'boundingbox'), offsetxy)
        }
    }

    #################################################
    ## standardize user model and parameterisation
    #################################################

    if ('formula' %in% class(model)) model <- list(model)
    model <- stdform (model)  ## named, no LHS
    if (CL) model$D <- NULL
    if (detector(traps(capthist)) %in% 'telemetry') model$g0 <- NULL
    details$param <- new.param(details, model, CL)
    ## intercept and fix certain models with bad defaults
    model <- updatemodel(model, detectfn, 9, c('g0', 'sigma'), c('b0', 'b1'))
    model <- updatemodel(model, detectfn, 10:13, c('g0', 'sigma'), c('beta0','beta1'))
    model <- updatemodel(model, detectfn, 14:18, 'g0', 'lambda0')

    #################################################
    ## which real parameters are fixed?
    #################################################

    ## c fixed by default in sigmak parameterisation
    if (details$param %in% 4:6) {
        if (! ("c" %in% names(model))) {
            ## default to fixed c = 0
            if (!("c" %in% names(fixed)))
                fixed$c <- 0
        }
    }
    fnames <- names(fixed)

    #################################################
    ## build default model and update with user input
    #################################################

    defaultmodel <- list(D=~1, g0=~1, lambda0=~1,  esa=~1, a0=~1,
                          sigma=~1, sigmak=~1, z=~1, w=~1, c=~1,
                          noneuc=~1, beta0=~1, beta1=~1,
                          sdS=~1, b0=~1, b1=~1, pID=~1, pmix=~1)
    defaultmodel <- replace (defaultmodel, names(model), model)

    #################################################
    # finite mixtures - 2009 12 10, ... 2014-10-16
    #################################################

    nmix <- get.nmix(model, capthist, hcov)
    if (nmix > 3)
        stop ("number of latent classes exceeds 3")
    if ((nmix>1) & !is.null(hcov) & !is.null(groups))
        stop ("hcov mixture model incompatible with groups")
    if ((nmix == 1) & ('pmix' %in% c(fnames,names(model))))
        stop ("pmix specified for invariant detection model")

    if ((nmix>1) & !('pmix' %in% fnames)) {
      if (is.null(model$pmix)) model$pmix <- ~1
      pmixvars <- all.vars(model$pmix)
      if (!any (pmixvars %in% c('h2','h3'))) ## add mixing h2 or h3
      {
        defaultmodel$pmix <- if (nmix == 2)
          update(model$pmix, ~. + h2)
        else
          update(model$pmix, ~. + h3)
      }
      else {
        defaultmodel$pmix <- model$pmix   ## use as-is
        badvar <- !(pmixvars %in% c('session','Session',sessioncov,'h2','h3'))
        if (any(badvar))
          stop ("formula for pmix may not include ", pmixvars[badvar])
      }
    }
    details$nmix <- nmix

    #################################################
    ## parameter names
    #################################################

    pnames <- valid.pnames (details, CL, detectfn, alltelem, sighting, nmix)

    #################################################
    ## test for irrelevant parameters in user's model
    #################################################

    OK <- names(model) %in% pnames
    if (any(!OK))
      stop ("parameters in model not consistent with detectfn etc. : ",
            paste(names(model)[!OK], collapse = ', '))
    OK <- fnames %in% pnames
    if (any(!OK))
      stop ("attempt to fix parameters not in model : ",
            paste(fnames[!OK], collapse = ', '))

    #################################################
    ## finalise model
    #################################################

    pnames <- pnames[!(pnames %in% fnames)]   ## drop fixed real parameters
    model <- defaultmodel[pnames]             ## select real parameters
    valid.model(model, CL, detectfn, hcov, details$userdist, names(sessioncov))
    vars <-  unlist(lapply(model, all.vars))

    #################################################
    ## Specialisations
    #################################################
    if (CL & !is.null(groups)) {
      groups <- NULL
      warning ("groups not valid with CL; groups ignored")
    }
    if (CL && var.in.model('g', model))
      stop ("'g' is not a valid effect when 'CL = TRUE'")
    if ((length(model) == 0) & !is.null(fixed))
      stop ("all parameters fixed")     ## assume want only LL

    ## mark-resight
    if ('pID' %in% names(model)) {
        pIDvars <- all.vars(model[['pID']])
        if (!all(pIDvars %in% c('session','Session',
                            names(sessioncov),
                            names(covariates(traps)))
                 ))
            warning ("predictors in model for pID may be invalid")
    }
    #################################################
    # Link functions (model-specific)
    #################################################

    defaultlink <- list(D='log', g0='logit', lambda0='log', esa='log',
                        a0='log', sigma='log', sigmak='log', z='log',
                        w='log', c='identity', noneuc='log',
                        beta0='identity', beta1='neglog', sdS='log',
                        b0='log', b1='neglog',  pID='logit',
                        pmix='logit', cut='identity')

    if (anycount) defaultlink$g0 <- 'log'
    link <- replace (defaultlink, names(link), link)
    link[!(names(link) %in% c(fnames,pnames))] <- NULL

    ##############################################
    # Prepare detection design matrices and lookup
    ##############################################
    memo ('Preparing detection design matrices', details$trace)
    design <- secr.design.MS (capthist, model, timecov, sessioncov, groups, hcov,
                              dframe, ignoreusage = details$ignoreusage)
    design0 <- secr.design.MS (capthist, model, timecov, sessioncov, groups, hcov,
                               dframe, ignoreusage = details$ignoreusage, naive = T,
                               bygroup = !CL)

    ############################################
    # Prepare density design matrix
    ############################################
    D.modelled <- !CL & is.null(fixed$D)
    smoothsetup <- list(D = NULL, noneuc = NULL)
    if (!D.modelled) {
       designD <- matrix(nrow = 0, ncol = 0)
       grouplevels <- 1    ## was NULL
       attr(designD, 'dimD') <- NA
       nDensityParameters <- integer(0)
    }
    else {
        grouplevels  <- group.levels(capthist,groups)
        if (!is.null(details$userDfn)) {
            ## may provide a function used by getD in functions.R
            ## userDfn(mask, beta[parindx$D], ngrp, nsession)
            designD <- details$userDfn
            if (!is.function(designD))
                stop ("details$userDfn should be a function")
            ## this form of call returns only coefficient names
            Dnames <- designD('parameters', mask)
        }
        else {
            memo ('Preparing density design matrix', details$trace)
            if (!all (all.vars(model$D) %in%
                      c('session', 'Session','g')) & details$param %in% c(4,5)) {
                ## 2014-09-30, 2014-10-04
                if (is.null(details$userdist))
                stop ("only session and group models allowed for density when details$param = ",
                      details$param)
            }
            temp <- D.designdata( mask, model$D, grouplevels, sessionlevels, sessioncov)
            if (any(smooths(model$D)))
                smoothsetup$D <- gamsetup(model$D, temp)
            ## otherwise, smoothsetup$D remains NULL
            designD <- general.model.matrix(model$D, temp, smoothsetup$D)
            attr(designD, 'dimD') <- attr(temp, 'dimD')

            Dnames <- colnames(designD)
        }
        nDensityParameters <- length(Dnames)
    }

    ############################################
    # Prepare non-Euclidean design matrix
    ############################################
    ## NE.modelled <- is.function(details$userdist)
    NE.modelled <- ('noneuc' %in% getuserdistnames(details$userdist)) &
        is.null(fixed$noneuc)
    if (!NE.modelled) {
       designNE <- matrix(nrow = 0, ncol = 0)
       grouplevelsNE <- 1    ## was NULL
       attr(designNE, 'dimD') <- NA
       nNEParameters <- integer(0)
    }
    else {
        grouplevelsNE  <- group.levels(capthist,groups)
        memo ('Preparing non-Euclidean parameter design matrix', details$trace)
        temp <- D.designdata( mask, model$noneuc, grouplevelsNE, sessionlevels, sessioncov)
        if (any(smooths(model$noneuc)))
            smoothsetup$noneuc <- gamsetup(model$noneuc, temp)
        ## otherwise, smoothsetup$NE remains NULL
        designNE <- general.model.matrix(model$noneuc, temp, smoothsetup$noneuc)
        attr(designNE, 'dimD') <- attr(temp, 'dimD')
        NEnames <- colnames(designNE)
        nNEParameters <- length(NEnames)
    }

    ############################################
    # Parameter mapping (general)
    ############################################
    if (!is.null(design$designMatrices))
        np <- sapply(design$designMatrices, ncol)
    else {
        np <- c(detectpar = 0)
    }
    np <- c(D = nDensityParameters, np, noneuc = nNEParameters)
    NP <- sum(np)
    parindx <- split(1:NP, rep(1:length(np), np))
    names(parindx) <- names(np)[np>0]
    if (!D.modelled) parindx$D <- NULL
    if (!NE.modelled) parindx$noneuc <- NULL

    ############################################
    # Optionally start from previous fit
    ############################################

    if (inherits(start, 'secr')) {
        ## use 'mapbeta' from score.test.R
        oldbeta <- coef(start)$beta
        if (details$nsim > 0)
            start <- oldbeta    ## chat simulations
        else {
            names(oldbeta) <- start$betanames
            oldnam <- start$betanames
            start <- mapbeta(start$parindx, parindx, oldbeta, NULL)
            if (!is.null(details$miscparm)) {
                nb <- length(start)
                start <- c(start, details$miscparm)
                oldnam <- oldnam[oldnam %in% names(details$miscparm)]
                start[oldnam] <- oldbeta[oldnam]
            }
        }
    }

    ############################################
    # send data to worker processes
    # do it once, not each evaluation
    ############################################
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = FALSE, useXDR = .Platform$endian=='big')
        ## clusterEvalQ(clust, requireNamespace('secr'))
         clusterExport(clust, c("capthist", "mask", "groups",
             "design","design0", "detectfn"), envir = environment())
    }
    else
        clust <- NULL

    ############################################
    # Single evaluation option
    ############################################
    .localstuff$iter <- 0
    if (details$LLonly) {
        if (is.null(start) | is.list(start))
            stop ("must provide vector of beta parameter values in 'start'")
        LL <- - secr.loglikfn (beta = start,
                               parindx    = parindx,
                               link       = link,
                               fixedpar   = fixed,
                               designD    = designD,
                               designNE   = designNE,
                               design     = design,
                               design0    = design0,
                               capthist   = capthist,
                               mask       = mask,
                               detectfn   = detectfn,
                               CL         = CL,
                               hcov       = hcov,
                               groups     = groups,
                               details    = details,
                               logmult    = TRUE,     ## add if possible
                               ncores     = ncores,
                               clust      = clust
        )

        return(c(logLik=LL))
    }
    ############################################

    ## estimate overdispersion by simulation
    ## if (nsim > 0) and details$chatonly then exit, returning only chat
    if (is.null(details$chat))
        details$chat <- matrix(1, nrow = length(sessionlevels), ncol = 2)
    else
        details$chat <- rep(details$chat,2)[1:2]  ## duplicate scalar
    if (details$nsim > 0) {
        TuTm <- function(x) !(is.null(Tu(x)) & is.null(Tm(x)))
        anysightings <- if (MS)
            any (sapply(capthist, TuTm))
        else TuTm(capthist)
        if (anysightings) {
            memo('Simulating sightings to estimate overdispersion...', details$trace)

            chat <- secr.loglikfn (beta       = start,
                                   parindx    = parindx,
                                   link       = link,
                                   fixedpar   = fixed,
                                   designD    = designD,
                                   designNE   = designNE,
                                   design     = design,
                                   design0    = design0,
                                   capthist   = capthist,
                                   mask       = mask,
                                   detectfn   = detectfn,
                                   CL         = CL,
                                   hcov       = hcov,
                                   groups     = groups,
                                   details    = details,
                                   logmult    = FALSE,
                                   ncores     = ncores,
                                   clust      = clust
            )
            if (details$chatonly)
                return(chat)   ## and no more!
            else {
                details$chat <- chat
                details$nsim <- 0
                ## and proceed to call secr.loglikfn again
            }
        }
    }

    ############################################
    # Start values (model-specific)
    # 'start' is vector of beta values (i.e. transformed)
    # or a list (secr >= 2.9.1)
    ############################################
    ## 2014-10-14
    ## if (is.null(start)) {
    ## 2015-01-06 only use autoini if required
    if (is.null(start) | is.list(start)) {
        start3 <- list(D = NA, g0 = NA, sigma = NA)
        ch <- if (MS) capthist[[details$autoini]] else capthist
        msk <- if (MS) mask[[details$autoini]] else mask

        requireautoini <- is.null(start) | !all(names(parindx) %in% names(start))
        if (requireautoini) {
            ## not for signal attenuation
            if (!(detectfn %in% c(9,10,11,12,13)) && !anypoly && !anytrans) {
                memo('Finding initial parameter values...', details$trace)
                ## autoini uses default buffer dbar * 4
                ## 2015-01-06
                if (nrow(ch)<5)
                    stop ("too few values in session 1 to determine start; set manually")
                start3 <- autoini (ch, msk, binomN = details$binomN,
                                   ignoreusage = details$ignoreusage)

                if (any(is.na(unlist(start3)))) {
                    warning ("'secr.fit' failed because initial values not found",
                             " (data sparse?); specify transformed values in 'start'")
                    return (list(call = cl, fit = NULL))
                }
                if (details$unmash & !CL) {
                    nmash <- attr(ch, 'n.mash')
                    if (!is.null(nmash)) {
                        n.clust <- length(nmash)
                        start3$D <- start3$D / n.clust
                    }
                }
                nms <- c('D', 'g0', 'sigma')
                nms <- paste(nms, '=', round(unlist(start3),5))
                memo(paste('Initial values ', paste(nms, collapse=', ')),
                     details$trace)
            }
            else warning ("using default starting values")
        }
        #--------------------------------------------------------------
        # assemble start vector
        ## revised 2014-12-04 to avoid sessions with no detections
        rpsv <- try(RPSV(ch, TRUE), silent = TRUE)
        if (inherits(rpsv, 'try-error'))
            rpsv <- NA
        n <- nrow(ch)

        default <- list(
            D       = ifelse (is.na(start3$D), 1, start3$D),
            g0      = ifelse (is.na(start3$g0), 0.1, start3$g0),
            lambda0 = -log(1-ifelse (is.na(start3$g0), 0.1, start3$g0)),
            sigma   = ifelse (is.na(start3$sigma), rpsv, start3$sigma),
            z       = 5,
            w       = 10,
            pID     = 0.7,
            noneuc  = 50,
            beta0   = details$cutval + 30,
            beta1   = -0.2,
            sdS     = 2,
            b0      = 2,      ## changed from 15 2010-11-01
            b1      = -0.1,
            pmix    = 0,      ## superceded below
            esa     = ifelse (is.na(start3$D), n / 5, n / start3$D),
            a0      = ifelse (is.na(start3$g0), 0.1 * rpsv^2, start3$g0 *
                      start3$sigma^2) / 10000 * 2 * pi
        )

        if (details$param %in% 4:6) {
            default$sigmak <- default$sigma * default$D^0.5
            default$c <- 0 ## but problems if take log(c)
        }

        if (detectfn %in% c(6)) {
            default$w <- default$sigma
            default$sigma <- default$sigma/2
        }
        if (detectfn %in% c(7)) {
            default$z <- default$sigma/5
        }
        if (detectfn %in% c(8, 18)) {
            default$z <- 1    ## cumulative gamma
        }
        if (anypoly | anytrans) {

            ## 2015-05-10 streamlined
            default$D <- 2 * nrow(ch) / maskarea(msk)

            ## 2015-10-11 better lambda0 -> g0 than vv...
            default$lambda0 <- sum(ch) / nrow(ch) / ncol(ch)
            
            ## default$g0 <- 1 - exp(-default$lambda)
            ## bugfix 2015-12-09
            default$g0 <- 1 - exp(-default$lambda0)
            
            if (details$binomN > 1)
                default$g0 <- default$g0 / details$binomN
            if ((details$binomN == 1) & (detector(traps(ch)) %in%
                                             c('polygon','transect'))) {
                ## assume using usage for binomN
                usge <- usage(traps(ch))
                default$g0 <- default$g0 / mean(usge[usge>0])
            }
            default$sigma <- rpsv
        }

        if (is.na(default$sigma)) default$sigma <- 20
        getdefault <- function (par) {
            transform (default[[par]], link[[par]])
        }

        ## 2014-10-14, 2014-11-12
        if (is.list(start)) {
            startnames <- names(start)
            default <- replace(default, startnames, start)
        }
        else startnames <- NULL

        start <- rep(0, NP)
        for ( i in 1:length(parindx) )
            start[parindx[[i]][1]] <- getdefault (names(model)[i])
        if ((details$nmix>1) & !('pmix' %in% fnames) & !('pmix' %in% startnames))
            start[parindx[['pmix']][1]] <- clean.mlogit((1:nmix)-0.5)[2]

        if (detectfn %in% c(12,13))
            start <- c(start, 46,3)    ## muN, sdN

        # D/ngrp when figure out where to calculate this

        if (sighting & is.null(fixed$pID))
             start[parindx$pID[1]] <- getdefault('pID')

        # start vector completed
        #--------------------------------------------------------------
    }
    ############################################
    ## ad hoc fix for experimental parameters
    ############################################
    if (is.null(details$miscparm))
        nmiscparm <- 0
    else {
        nmiscparm <- length(details$miscparm)
        if (length(start) < max(unlist(parindx)) + nmiscparm )
            start <- c(start, details$miscparm)
    }
    if (detector(traps(capthist)) %in% c('signalnoise'))
        nmiscparm <- 2

    NP <- NP + nmiscparm
    stopifnot (length(start) == NP)

    ############################################
    # Fixed beta parameters
    ############################################
    fb <- details$fixedbeta
    if (!is.null(fb)) {
        if (!(length(fb)== NP))
            stop ("invalid fixed beta - require NP-vector")
        if (sum(is.na(fb))==0)
            stop ("cannot fix all beta parameters")
        ## drop unwanted betas; remember later to adjust parameter count
        start <- start[is.na(fb)]
        NP <- length(start)
    }
    ############################################
    # Variable names (general)
    ############################################
    betanames <- unlist(sapply(design$designMatrices, colnames))
    names(betanames) <- NULL
    realnames <- names(model)
    if (D.modelled) betanames <- c(paste('D', Dnames, sep='.'), betanames)
    if (NE.modelled) betanames <- c(betanames, paste('noneuc', NEnames, sep='.'))
    betanames <- sub('..(Intercept))','',betanames)

    ############################################
    # Variable names (model-specific)
    ############################################

    if (detectfn %in% c(12,13)) {
        betanames <- c(betanames, 'muN', 'sdN')
        realnames <- c(realnames, 'muN', 'sdN')
    }
    else if (nmiscparm>0) {
        miscnames <- names(details$miscparm)
        if (is.null(miscnames))
            miscnames <- paste('miscparm', 1:nmiscparm, sep='')
        betanames <- c(betanames, miscnames)
    }

    ## allow for fixed beta parameters
    if (!is.null(details$fixedbeta))
        betanames <- betanames[is.na(details$fixedbeta)]
    betaw <- max(max(nchar(betanames)),8)   # for 'trace' formatting

    ############################################
    # Maximize likelihood
    ############################################

    lcmethod <- tolower(method)
    if (lcmethod != 'none') {
        memo('Maximizing likelihood...', details$trace)
        if (details$trace)
            cat('Eval     Loglik', formatC(betanames, format='f', width=betaw), '\n')
    }
    loglikefn <- secr.loglikfn

    ## arguments always passed to loglikefn
    secrargs <- list(
                     link       = link,
                     fixedpar   = fixed,
                     parindx    = parindx,
                     capthist   = capthist,
                     mask       = mask,
                     CL         = CL,
                     detectfn   = detectfn,
                     designD    = designD,
                     designNE   = designNE,
                     design     = design,
                     design0    = design0,
                     hcov       = hcov,
                     groups     = groups,
                     details    = details,
                     logmult    = TRUE,     # add if possible
                     ncores     = ncores,
                     clust      = clust,
                     betaw      = betaw)    # for trace format
    ############################################
    ## calls for specific maximisation methods
    ## 2013-04-21
    if (NP == 1) {
        lcmethod <- "optimise"
        signs <- c(-1,1) * sign(start)
        args <- list (f         = loglikefn,
                      interval  = start * (1 + details$intwidth2 * signs))
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- try(do.call (optimise, args))
        if (inherits(this.fit, 'try-error'))
            warning ("univariate search for minimum failed")
        this.fit$par <- this.fit$minimum
        this.fit$value <- this.fit$objective
        if (details$hessian != "none")
            details$hessian <- "fdHess"
    }
    else if (lcmethod %in% c('newton-raphson')) {

        args <- list (p         = start,
                      f         = loglikefn,
                      hessian   = tolower(details$hessian)=='auto',
                      stepmax   = 10)
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))

        this.fit <- do.call (nlm, args)

        this.fit$par <- this.fit$estimate     # copy for uniformity
        this.fit$value <- this.fit$minimum    # copy for uniformity
        if (this.fit$code > 2)
            warning ("possible maximization error: nlm returned code ",
                this.fit$code, ". See ?nlm")
    }
    #-----------------------------------------------------------------
    else if (method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B",
                             "SANN", "Brent")) {
        args <- list(par     = start,
                     fn      = loglikefn,
                     hessian = tolower(details$hessian)=='auto',
                     method  = method)
        args <- c(args, secrargs)
        args <- replace (args, names(list(...)), list(...))
        this.fit <- do.call (optim, args)
        if (this.fit$convergence != 0)
            warning ("probable maximization error: optim returned convergence ",
                this.fit$convergence, ". See ?optim")
    }
    #-----------------------------------------------------------------
    # Hessian-only 2013-02-23
    else if (lcmethod %in% 'none') {
        memo ('Computing Hessian with fdHess in nlme', details$trace)
        loglikfn <- function (beta) {
            do.call(secr.loglikfn, c(list(beta=beta), secrargs))
        }
        grad.Hess <- nlme::fdHess(start, fun = loglikfn, .relStep = 0.001, minAbsPar=0.1)
        this.fit <- list (value = loglikfn(start), par = start,
                          gradient = grad.Hess$gradient,
                          hessian = grad.Hess$Hessian)
        biasLimit <- NA   ## no bias assessment
    }
    else stop ("maximization method", method, "not recognised")
    ############################################################################

    this.fit$method <- method         ## remember which method we used...
    covar <- NULL
    N <- NULL
    if (this.fit$value > 1e9) {     ## failed
        # this.fit$beta[] <- NA
        this.fit$par[] <- NA
    }
    else {

    ############################################
    ## Variance-covariance matrix
    ############################################

        if (tolower(details$hessian)=='fdhess') {
                memo ('Computing Hessian with fdHess in nlme', details$trace)
                loglikfn <- function (beta) {
                    do.call(secr.loglikfn, c(list(beta=beta), secrargs))
                }
                grad.Hess <- nlme::fdHess(this.fit$par, fun = loglikfn,
                                          .relStep = 0.001, minAbsPar = 0.1)
                this.fit$hessian <- grad.Hess$Hessian
        }
        if (!is.null(this.fit$hessian)) {
            covar <- try(MASS::ginv(this.fit$hessian))
            if (inherits(covar, "try-error")) {
                warning ("could not invert Hessian to compute ",
                         "variance-covariance matrix")
                covar <- matrix(nrow = NP, ncol = NP)
            }
            else if (any(diag(covar)<0)) {
                warning ("at least one variance calculation failed ")
            }
            dimnames(covar) <- list(betanames, betanames)
        }

        ## predicted D across mask
        if (!CL) {
            D <- getD (designD, this.fit$par, mask, parindx, link, fixed,
                       grouplevels, sessionlevels, 'D')
            N <- t(apply(D, 2:3, sum, drop = FALSE))

            if (inherits(mask, 'linearmask')) {
                if (ms(mask)) {
                    scale <- sapply(mask, attr, 'spacing')
                }
                else {
                    scale <- attr(mask, 'spacing')
                }
                scale <- scale / 1000   ## per km
            }
            else {
                if (ms(mask)) {
                    scale <- sapply(mask, attr, 'area')
                }
                else {
                    scale <- attr(mask, 'area')
                }
            }
            ## rows are sessions
            N <- sweep(N, FUN = '*', MARGIN = 1, STATS = scale)
        }
    }

    ############################################
    ## form output list
    ############################################

    desc <- packageDescription("secr")  ## for version number

    ## if density modelled then smoothsetup already contains D,noneuc
    ## otherwise an empty list; now add detection parameters from
    ## secr.design.MS
    smoothsetup <- c(smoothsetup, design$smoothsetup)

    output <- list (call = cl,
                  capthist = capthist,
                  mask = mask,
                  detectfn = detectfn,
                  CL = CL,
                  timecov = timecov,
                  sessioncov = sessioncov,
                  hcov = hcov,
                  groups = groups,
                  dframe = dframe,
                  design = design,
                  design0 = design0,
                  start = start,
                  link = link,
                  fixed = fixed,
                  parindx = parindx,
                  model = model,
                  details = details,
                  vars = vars,
                  betanames = betanames,
                  realnames = realnames,
                  fit = this.fit,
                  beta.vcv = covar,
                  smoothsetup = smoothsetup,
                  N = N,
                  version = desc$Version,
                  starttime = starttime,
                  proctime = (proc.time() - ptm)[3]
             )

    class(output) <- 'secr'

    if (usebuffer & !is.na(biasLimit)) {
        test <- try(bufferbiascheck(output, buffer, biasLimit))
    }

    memo(paste('Completed in ', round(output$proctime,2), ' seconds at ',
        format(Sys.time(), "%H:%M:%S %d %b %Y"),
        sep=''), details$trace)

    if (ncores > 1) {
        stopCluster(clust)
    }

    output

}
################################################################################
