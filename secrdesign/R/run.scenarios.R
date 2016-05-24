#c#############################################################################
## package 'secrdesign'
## run.scenarios.R
## 2013-02-23, 24, 25, 26, 28
## 2014-02-07, 2014-02-09, 2014-02-10
## currently only for single-session model
## 2014-04-06 drop logfile argument, reorder arguments, add pop.args
## 2014-04-14 fit.models
## 2014-04-26 IHP
## 2014-04-27 automatically wrap mask, pop.args, det.args if not in list form
## 2014-09-03 linear mask tweaks
## 2014-11-23 groups
## 2014-11-24 new functions makeCH and processCH used by onesim()
## 2015-01-26 defaultextractfn updated
## 2015-01-26 streamlined outputtype (function getoutputtype called by both run.scenarios and fit.models)
## 2015-01-26 scen and repl argument for fit.models
## 2015-01-27 more robust handling of start values in processCH
## 2015-11-03 default extract fn reports unmarked and nonID sightings
## 2015-11-03 adapted for sim.resight etc.

###############################################################################
wrapifneeded <- function (args, default) {
    if (any(names(args) %in% names(default)))
        list(args)  ## assume single; wrap
    else
        args
}
###############################################################################
## complete a partial argument list (arg) with default args
## always return a list of arg lists long enough to match max(index)
fullargs <- function (args, default, index) {
    if (is.null(args)) {
        full.args <- list(default)
    }
    else {
        nind <- max(index)
        if (length(args) < nind) stop("too few components in args")
        tmpargs <- vector('list', nind)
        for (i in 1:nind) {
            if (is.character(args[[i]]))
                args[[i]] <- match.arg(args[[i]], default[[i]])
            ## naked fn gives trouble here... 2014-09-03
            tmpargs[[i]] <- replace (default, names(args[[i]]), args[[i]])
            if (is.character(tmpargs[[i]]))
                tmpargs[[i]] <- tmpargs[[i]][1]
        }
        full.args <- tmpargs
    }
    full.args
}
###############################################################################

defaultextractfn <- function(x) {
    ## 2015-01-27
    if (inherits(x, 'try-error')) {
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
    }
    else if (inherits(x, 'capthist')) {
        ## summarised raw data
        counts <- function(CH) {
            ## for single-session CH
            if (nrow(CH)==0)  ## 2015-01-24
                c(n=0, ndet=0, nmov=0, dpa = NA)
            else {
                nmoves <- sum(unlist(sapply(moves(CH), function(y) y>0)))
                ## detectors per animal
                dpa <- if (length(dim(CH)) == 2)
                    mean(apply(abs(CH), 1, function(y) length(unique(y[y>0]))))
                else
                    mean(apply(apply(abs(CH), c(1,3), sum)>0, 1, sum))
                if (sighting(traps(CH))) {
                    unmarked <- if (is.null(Tu <- Tu(CH))) NA else sum(Tu)
                    nonID <- if (is.null(Tm <- Tm(CH))) NA else sum(Tm)
                    nzero <- sum(apply(abs(CH),1,sum) == 0)
                    c(n=nrow(CH), ndet=sum(abs(CH)>0), nmov=nmoves, dpa = dpa,
                      unmarked=unmarked, nonID = nonID, nzero = nzero)
                }
                else
                    c(n=nrow(CH), ndet=sum(abs(CH)>0), nmov=nmoves, dpa = dpa)
            }
        }
        if (ms(x))
            unlist(lapply(x, counts))
        else {
            gp <- covariates(x)$group
            if (is.null(gp))
                counts(x)
            else
                unlist(lapply(split(x,gp,dropnullocc=TRUE), counts))
        }
    }
    else if (inherits(x,'secr') & (!is.null(x$fit))) {
        ## fitted model:
        ## default predictions of 'real' parameters
        out <- predict(x)
        if (is.data.frame(out))
            out
        else {
            ## 2015-01-26
            warning ("summarising only first session, group or mixture class")
            out[[1]]
        }
    }
    else
        ## null output: dataframe of 0 rows and 0 columns
        data.frame()
}

#####################
makeCH <- function (scenario, trapset, full.pop.args, full.det.args, mask, multisession) {
    ns <- nrow(scenario)
    with( scenario, {
        CH <- vector(mode = 'list', ns)

        for (i in 1:ns) {
            #####################
            ## retrieve data
            grid <- trapset[[trapsindex[i]]]
            poparg <- full.pop.args[[popindex[i]]]
            detarg <- full.det.args[[detindex[i]]]

            #####################
            ## override D, core, buffer
            if (inherits(mask, 'linearmask'))               ## force to linear...
                poparg$model2D <- 'linear'
            if (poparg$model2D %in% c('IHP', 'linear')) {   ## linear
                poparg$core <- mask
                ## for 'linear' case we may want a constant numeric value
                if (!is.character(poparg$D) & (length(poparg$D)<nrow(mask)))
                    poparg$D <- D[i]
                if (nrepeats[i]!=1)
                    stop("nrepeats > 1 not allowed for IHP, linear")
            }
            else {
                poparg$core <- attr(mask, "boundingbox")
                poparg$D <- D[i]*nrepeats[i]
            }
            poparg$buffer <- 0

            #####################
            ## generate population
            pop <- do.call(sim.popn, poparg)

            #####################
            ## form dp for sim.capthist or sim.resight
            ## form par for starting values in secr.fit()
            ## 'par' does not allow for varying link or any non-null model (b, T etc.)
            if (detectfn[i] %in% 14:18) {
                dp <- list(lambda0 = lambda0[i], sigma = sigma[i],
                           recapfactor = recapfactor[i])
            }
            else {
                dp <- list(g0 = g0[i], sigma = sigma[i],
                           recapfactor = recapfactor[i])
            }

            ## 2016-03-06
            dp <- replace (dp, names(detarg$detectpar), detarg$detectpar)

            #####################
            ## override det args as required
            detarg$traps <- grid
            detarg$popn <- pop
            detarg$detectpar <- dp
            detarg$detectfn <- detectfn[i]
            if (!is.null(markocc(grid))) {
                detarg$noccasions <- length(markocc(grid))
                if (detarg$noccasions != noccasions[i])
                    warning("length of markocc attribute overrides noccasions in scenario")
            }
            else
                detarg$noccasions <- noccasions[i]

            #####################
            ## simulate detection
            if (sighting(grid)) {
                CHi <- do.call(sim.resight, detarg)
            }
            else
                CHi <- do.call(sim.capthist, detarg)
            if (!is.na(nrepeats[i]))
                attr(CHi, "n.mash") <- rep(NA, nrepeats[i])
            CH[[i]] <- CHi
        }
        if (ns > 1) {
            ## assume a 'group' column is present if ns>1
            names(CH) <- 1:ns
            if (multisession) {
                CH <- MS.capthist(CH)
                if (!is.null(group)) session(CH) <- group
                CH
            }
            else {
                nc <- sapply(CH, nrow)
                CH$verify <- FALSE
                CH <- do.call(rbind.capthist, CH)
                covariates(CH)$group <- rep(group, nc)
                CH
            }
        }
        else {
            CH[[1]]
        }
    })
}
#####################
processCH <- function (scenario, CH, fitarg, extractfn, fit, ...) {
    if (!fit) {
        extractfn(CH, ...)
    }
    else {
        ## form par for starting values in secr.fit()
        ## 'par' does not allow for varying link or any non-null model (b, T etc.)
        ## D, lambda0, g0, sigma are columns in 'scenario'

        par <- with(scenario, {
            wt <- D/sum(D)
            ## 2015-01-27
            if (detectfn[1] %in% 14:18) {
                list(D = sum(D), lambda0 = sum(lambda0*wt), sigma = sum(sigma*wt))
            }
            else {
                list(D = sum(D), g0 = sum(g0*wt), sigma = sum(sigma*wt))
            }
        })
        ## prepare arguments for secr.fit()
        fitarg$capthist <- CH

        if (is.null(fitarg$model))
            fitarg$model <- defaultmodel(fitarg$CL, fitarg$detectfn)
        if (fitarg$start[1] == 'true') {
            ## check to see if simple 'true' values will work
            ## requires intercept-only model for all parameters
            model <- eval(fitarg$model)
            if (!is.list(model)) model <- list(model)
            vars <- unlist(lapply(lapply(model, terms), attr, 'term.labels'))
            if (fitarg$CL) par$D <- NULL
            if ('h2' %in% vars) par$pmix <- 0.5
            fitarg$start <- par
            if ((length(vars) != 0) & (fitarg$method == 'none')) {
                ## not yet ready for interspersed beta coef
                stop("method = 'none' requires full start vector for complex models")
            }
        }
        fitarg$trace <- FALSE

        ##-------------------------------------------------------------------
        ## 2015-11-03 code for overdispersion adjustment of mark-resight data
        ## no simulations in first iteration; defer hessian
        chatnsim <- fitarg$details$nsim
        if (abs(chatnsim)>0) {
            fitarg$details <- as.list(replace(fitarg$details, 'nsim', 0))
            fitarg$details <- as.list(replace(fitarg$details, 'hessian', FALSE))
        }
        ##-------------------------------------------------------------------

        fit <- try(do.call(secr.fit, fitarg))

        ##-------------------------------------------------------------------
        ## 2015-11-03 code for overdispersion adjustment of mark-resight data
        ## 2016-03-06 patched
        if (sighting(traps(CH)) & !inherits(fit, 'try-error')) {
            if ((abs(chatnsim) > 0) &  (logLik(fit)>-1e9)) {
                fitarg$details$nsim <- abs(chatnsim)
                fitarg$details$hessian <- TRUE
                fit$call <- NULL
                fitarg$start <- fit
                if (chatnsim<0)
                    fitarg$method <- "none"
                fit <- try(do.call(secr.fit, fitarg))
            }
        }
        ##-------------------------------------------------------------------

        extractfn(fit, ...)
    }
}
#####################

getoutputtype <- function (output) {
    typical <- output[[1]][[1]]
    ## outputtype: secrfit, predicted, coef, numeric
    outputtype <-
        if (inherits(typical, 'secr'))
            'secrfit'
        else if (inherits(typical, 'data.frame')) {
            if (all(c('estimate','SE.estimate','lcl','ucl') %in% names(typical)) &
                    any(c('R.N','E.N') %in% rownames(typical)))
                'regionN'
            else if ( all(c('estimate','SE.estimate','lcl','ucl', 'CVn') %in% names(typical)))
                'derived'
            else if (all(c('estimate','SE.estimate','lcl','ucl') %in% names(typical)))
                'predicted'
            else if (all(c('beta','SE.beta','lcl','ucl') %in% names(typical)))
                'coef'
            else 'user'
        }
        else if (inherits(typical, 'capthist'))          ## rawdata
            'capthist'
        else if (is.vector(typical, mode = 'numeric'))   ## usually, summary of unfitted data
            'selectedstatistics'
        else
            'user'
    outputtype
}
#####################

getoutputclass <- function (outputtype) {
    switch (outputtype,
            secrfit = c("fittedmodels", 'secrdesign', 'list'),
            predicted = c("estimatetables", 'secrdesign', 'list'),
            derived = c("estimatetables", 'secrdesign', 'list'),
            regionN = c("estimatetables", 'secrdesign', 'list'),
            coef = c("estimatetables", 'secrdesign', 'list'),
            user = c("estimatetables", 'secrdesign', 'list'),
            capthist = c("rawdata", 'secrdesign', 'list'),
            selectedstatistics = c("selectedstatistics", 'secrdesign', 'list'),
            list      ## otherwise
    )
}

###############################################################################
run.scenarios <- function (nrepl,  scenarios, trapset, maskset, xsigma = 4,
    nx = 32, pop.args, det.args, fit = FALSE, fit.args, chatnsim = 0, extractfn = NULL,
    multisession = FALSE, ncores = 1, seed = 123,  ...) {

    #--------------------------------------------------------------------------
    onesim <- function (scenario) {
        ## 2014-11-23 allow multi-line scenarios
        ## only one mask an fitarg allowed per scenario
        fitarg <- full.fit.args[[scenario$fitindex[1]]]
        fitarg$mask <- maskset[[scenario$maskindex[1]]]
        CH <- makeCH(scenario, trapset, full.pop.args, full.det.args,
                     fitarg$mask, multisession)
        processCH(scenario, CH, fitarg, extractfn, fit, ...)
    }
    #--------------------------------------------------------------------------
    runscenario <- function(x) {
        out <- vector('list', nrepl)
        for (r in 1:nrepl) out[[r]] <- onesim(x)
        message("Completed scenario ", x$scenario[1])
        flush.console()
        out
    }
    ##--------------------------------------------------------------------------

    ## mainline
    ## record start time etc.
    ptm  <- proc.time()
    cl   <- match.call(expand.dots = TRUE)
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")
    if (ncores > nrow(scenarios))
        stop ("ncores exceeds number of scenarios")
    if (multisession & !anyDuplicated(scenarios$scenario))
        warning ("multisession = TRUE ignored because no scenario duplicated")

    ##--------------------------------------------
    ## default extractfn
    if (is.null(extractfn)) {
        extractfn <- defaultextractfn
    }
    ##--------------------------------------------
    ## preprocess inputs
    if (inherits(trapset, 'traps'))   ## otherwise assume already list of traps
        trapset <- list(trapset)
    if (!missing(maskset)) {
        if (inherits(maskset, 'mask'))   ## otherwise assume already list of masks
            maskset <- list(maskset)
        if (is.null(names(maskset)))
            names(maskset) <- paste('mask',1:length(maskset), sep='')
    }
    nk <- length(trapset)
    if (is.null(names(trapset)))
        names(trapset) <- paste('traps',1:nk, sep='')

    dettype <- sapply(trapset, detector)[scenarios$trapsindex]

    ## mark-resight options 2015-11-03
    sight <- sapply(trapset, function(x)
        if(ms(x)) sighting(x[[1]]) else  sighting(x))

    if (!(all(sight) | all(!sight)))
        stop ("cannot mix sighting and nonsighting simulations")
    sight <- any(sight)

    OK <- !any((scenarios$nrepeats>1) & (dettype == "single"))
    OK <- if(is.na(OK)) TRUE else OK
    if (!OK)
        warning("single-catch traps violate independence assumption for nrepeats > 1")

    ##---------------------------------------------
    ## allow user changes to default sim.popn arguments
    default.args <- as.list(args(sim.popn))[1:12]
    default.args$model2D  <- eval(default.args$model2D)[1]   ## 2014-09-03
    if (missing(pop.args)) pop.args <- NULL
    pop.args <- wrapifneeded(pop.args, default.args)
    full.pop.args <- fullargs (pop.args, default.args, scenarios$popindex)

    ##---------------------------------------------
    ## allow user changes to default sim.capthist or sim.resight arguments
    if (sight)
        default.args <- as.list(args(sim.resight))[c(2,4:9)]  ## drop traps & ... argument
    else
        default.args <- as.list(args(sim.capthist))[1:15]

    if (missing(det.args)) det.args <- NULL
    det.args <- wrapifneeded(det.args, default.args)
    full.det.args <- fullargs (det.args, default.args, scenarios$detindex)

    ## load detect args back into scenarios?

    ##---------------------------------------------
    ## allow user changes to default secr.fit arguments
    default.args <- as.list(args(secr.fit))[1:21]
    default.args$biasLimit <- NA   ## never check
    default.args$verify <- FALSE   ## never check
    default.args$start <- "true"   ## known values
    default.args$detectfn <- 0     ## halfnormal
    if (missing(fit.args)) fit.args <- NULL
    fit.args <- wrapifneeded(fit.args, default.args)
    full.fit.args <- fullargs (fit.args, default.args, scenarios$fitindex)
    full.fit.args[[1]]$details <- as.list(replace(full.fit.args$details,'nsim',chatnsim))
    ##--------------------------------------------
    ## construct masks as required
    if (missing(maskset)) {
        trapsigma <- scenarios[, c('trapsindex', 'sigma'), drop = FALSE]
        uts <- unique (trapsigma)
        code <- apply(uts, 1, paste, collapse='.')
        scenarios$maskindex <- match(apply(trapsigma,1,paste,collapse='.'), code)
        maskset <- vector('list', nrow(uts))
        for (k in 1:nrow(uts))
            maskset[[k]] <- make.mask(trapset[[uts$trapsindex[k]]], buffer = uts$sigma[k] * xsigma,
                                   type = 'trapbuffer', nx = nx)
    }
    else {
        uts <- NULL
        if (is.null(scenarios$maskindex)) {
            if ((length(maskset) == length(trapset)))
                scenarios[,'maskindex'] <- scenarios$trapsindex
            else if (length(maskset) == 1)
                scenarios[,'maskindex'] <- 1
            else
                stop ("for irregular maskset provide maskindex as a column in scenarios")
        }
        if (max(scenarios$maskindex) > length(maskset))
            stop ("maskindex does not match maskset")
    }

    #--------------------------------------------
    ## override nrepeats and D in scenarios when IHP distribution
    for (i in 1:nrow(scenarios)) {
        pi <- scenarios$popindex[i]
        mi <- scenarios$maskindex[i]
        if ((full.pop.args[[pi]]$model2D %in% c('IHP', 'linear')) &
            (is.character(full.pop.args[[pi]]$D))) {          ## linear 2014-09-03
            avD <- mean (covariates(maskset[[mi]])[,full.pop.args[[pi]]$D])
            scenarios[i, 'nrepeats'] <- 1   ## override
            scenarios[i, 'D'] <- avD
        }
    }

    #--------------------------------------------
    ## run simulations
    tmpscenarios <- split(scenarios, scenarios$scenario)
    if (ncores > 1) {
        list(...)    ## ensures promises evaluated see parallel vignette 2015-02-02
        clust <- makeCluster(ncores, methods = TRUE)
        clusterSetRNGStream(clust, seed)
        output <- parLapply(clust, tmpscenarios, runscenario)
        stopCluster(clust)
    }
    else {
        set.seed (seed)
        output <- lapply(tmpscenarios, runscenario)
    }
    ##-------------------------------------------
    ## tidy output
    outputtype <- getoutputtype(output)
    if (outputtype == 'selectedstatistics')
        ## collapse replicates within a scenario into a matrix
        output <- lapply(output, do.call, what = rbind)
    message("Completed in ", round((proc.time() - ptm)[3]/60,3), " minutes")
    desc <- packageDescription("secrdesign")  ## for version number
    value <- list (call = cl,
                   version = paste('secrdesign', desc$Version),
                   starttime = starttime,
                   proctime = (proc.time() - ptm)[3],
                   scenarios = scenarios,
                   trapset = trapset,
                   maskset = if (is.null(uts)) maskset else NULL,
                   xsigma = xsigma,
                   nx = nx,
                   pop.args = pop.args,
                   det.args = det.args,
                   fit = fit,
                   fit.args = fit.args,
                   extractfn = extractfn,
                   seed = seed,
                   chatnsim = chatnsim,
                   nrepl = nrepl,
                   output = output,
                   outputtype = outputtype
    )
    class(value) <- getoutputclass(outputtype)
    if (outputtype == 'regionN')
        attr(value, 'regionsize') <- sapply(output, function(x) attr(x[[1]], 'regionsize'))

    value
}


########################################################################################

## version of run.scenarios that accepts existing data and
## expands scenarios for multiple model definitions

fit.models <- function (rawdata, fit = FALSE, fit.args, chatnsim = 0, extractfn = NULL,
                        ncores = 1, scen, repl, ...) {

    #--------------------------------------------------------------------------
    onesim <- function (scenario, CH) {
        fitarg <- full.fit.args[[scenario$fitindex[1]]]
        fitarg$mask <- maskset[[scenario$maskindex[1]]]
        processCH(scenario, CH, fitarg, extractfn, fit, ...)
    }
    #--------------------------------------------------------------------------
    runscenario <- function(x) {
        out <- vector('list', nrepl)
        for (r in 1:nrepl) {
            ## match by name, not number 2015-01-27
            scenID <- as.character(trunc(x$scenario[1]))
            out[[r]] <- onesim(x, CHlist[[scenID]][[r]])
        }
        message("Completed scenario ", x$scenario[1])
        flush.console()
        out
    }
    ##--------------------------------------------------------------------------
    ## mainline
    if (!inherits(rawdata, "rawdata"))
        stop ("requires rawdata output from run.scenarios()")

    ## optionally select which scenarios to fit
    if (missing(scen)) {
        scen <- unique(rawdata$scenarios$scenario)
    }
    else {
        scen <- unique(scen)
        if (!all(scen %in% unique(rawdata$scenarios$scenario)))
            stop ("invalid scen argument")
        if (length(scen)<1)
            stop ("invalid scen argument")
    }
    ## optionally select which replicates to fit
    if (missing(repl)) {
        nrepl <- rawdata$nrepl
        repl <- 1:nrepl
    }
    else {
        repl <- unique(repl)
        if (!all(repl %in% 1:rawdata$nrepl))
            stop ("invalid repl argument")
        nrepl <- length(repl)
        if (nrepl<1)
            stop ("invalid repl argument")
    }

    CHlist <- lapply(rawdata$output[scen], '[', repl)
    scenarios <- rawdata$scenarios[rawdata$scenarios$scenario %in% scen,]

    trapset <- rawdata$trapset
    maskset <- rawdata$maskset
    ## record start time etc.
    ptm  <- proc.time()
    cl   <- match.call(expand.dots = TRUE)
    starttime <- format(Sys.time(), "%H:%M:%S %d %b %Y")

    ##--------------------------------------------
    ## default extractfn
    if (is.null(extractfn)) {
        extractfn <- defaultextractfn
    }

    ##---------------------------------------------
    ## allow user changes to default secr.fit arguments
    default.args <- as.list(args(secr.fit))[1:21]
    default.args$biasLimit <- NA   ## never check
    default.args$verify <- FALSE   ## never check
    default.args$start <- "true"   ## known values
    default.args$detectfn <- 0     ## halfnormal
    if (missing(fit.args)) fit.args <- NULL
    fit.args <- wrapifneeded(fit.args, default.args)
    nfit <- length(fit.args)
    if (nfit > 1) {
        ## expand scenarios by the required number of different model fits
        ##      scenarios <- scenarios[rep(scenarios$scenario, each = nfit),]
        scenarios <- scenarios[rep(1:nrow(scenarios), each = nfit),]
        scenarios$fitindex <- rep(1:nfit, length.out = nrow(scenarios))
        ## assign new unique scenario number by adding decimal fraction
        scenarios$scenario <- scenarios$scenario + scenarios$fitindex /
            10 ^ trunc(log10(nfit)+1)
        scenarios <- scenarios[order(scenarios$scenario),]
        rownames(scenarios) <- 1:nrow(scenarios)
    }
    full.fit.args <- fullargs (fit.args, default.args, scenarios$fitindex)
    full.fit.args[[1]]$details <- as.list(replace(full.fit.args$details,'nsim',chatnsim))

    ## construct masks as required
    if (is.null(maskset)) {
        trapsigma <- scenarios[, c('trapsindex', 'sigma'), drop = FALSE]
        uts <- unique (trapsigma)
        code <- apply(uts, 1, paste, collapse='.')
        scenarios$maskindex <- match(apply(trapsigma,1,paste,collapse='.'), code)
        maskset <- vector('list', nrow(uts))
        for (k in 1:nrow(uts))
            maskset[[k]] <- make.mask(trapset[[uts$trapsindex[k]]], buffer = uts$sigma[k] *
                                      rawdata$xsigma, type = 'trapbuffer', nx = rawdata$nx)
    }
    else {
        uts <- NULL
        if (is.null(scenarios$maskindex)) {
            if ((length(maskset) == length(trapset)))
                scenarios[,'maskindex'] <- scenarios$trapsindex
            else if (length(maskset) == 1)
                scenarios[,'maskindex'] <- 1
            else
                stop ("for irregular maskset provide maskindex as a column in scenarios")
        }
        if (max(scenarios$maskindex) > length(maskset))
            stop ("maskindex does not match maskset")
    }

    if (ncores > nrow(scenarios))
        stop ("ncores exceeds number of scenarios")

    #--------------------------------------------
    ## run simulations
    tmpscenarios <- split(scenarios, scenarios$scenario)
    if (ncores > 1) {
        clust <- makeCluster(ncores, methods = TRUE)
        output <- parLapply(clust, tmpscenarios, runscenario)
        on.exit(stopCluster(clust))
    }
    else {
        output <- lapply(tmpscenarios, runscenario)
    }
    ##-------------------------------------------
    ## tidy output
    outputtype <- getoutputtype(output)
    if (outputtype == 'selectedstatistics')
        ## collapse replicates within a scenario into a matrix
        output <- lapply(output, do.call, what = rbind)
    message("Completed in ", round((proc.time() - ptm)[3]/60,3), " minutes")
    desc <- packageDescription("secrdesign")  ## for version number
    value <- list (call = cl,
                   version = paste('secrdesign', desc$Version),
                   starttime = starttime,
                   proctime = (proc.time() - ptm)[3],
                   scenarios = scenarios,
                   trapset = trapset,
                   maskset = maskset,
                   xsigma = rawdata$xsigma,
                   nx = rawdata$nx,
                   pop.args = rawdata$pop.args,
                   det.args = rawdata$det.args,
                   fit = fit,
                   fit.args = fit.args,
                   extractfn = extractfn,
                   seed = rawdata$seed,
                   nrepl = nrepl,     ## rawdata$nrepl, 2015-01-26
                   output = output,
                   outputtype = outputtype
                   )
    class(value) <- getoutputclass (outputtype)
    if (outputtype == 'regionN')
        attr(value, 'regionsize') <- sapply(output, function(x) attr(x[[1]], 'regionsize'))

    value
}
