################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Maximum Likelihood inference for the two-component spatio-temporal intensity
### model described in Meyer et al (2012), DOI: 10.1111/j.1541-0420.2011.01684.x
###
### Copyright (C) 2009-2016 Sebastian Meyer
### $Revision: 1617 $
### $Date: 2016-03-10 14:31:43 +0100 (Don, 10. Mär 2016) $
################################################################################


## model.frame() evaluates 'subset' and '...' with 'data'
utils::globalVariables(c("tile", "type", "BLOCK", ".obsInfLength", ".bdist",
                         "area"))

twinstim <- function (
    endemic, epidemic, siaf, tiaf, qmatrix = data$qmatrix,
    data, subset, t0 = data$stgrid$start[1], T = tail(data$stgrid$stop,1),
    na.action = na.fail, start = NULL, partial = FALSE,
    epilink = "log", control.siaf = list(F=list(), Deriv=list()),
    optim.args = list(), finetune = FALSE,
    model = FALSE, cumCIF = FALSE, cumCIF.pb = interactive(),
    cores = 1, verbose = TRUE
    )
{

    ####################
    ### Preparations ###
    ####################

    ptm <- proc.time()
    cl <- match.call()
    partial <- as.logical(partial)
    finetune <- if (partial) FALSE else as.logical(finetune)

    ## (inverse) link function for the epidemic linear predictor of event marks
    epilink <- match.arg(epilink, choices = c("log", "identity"))
    epilinkinv <- switch(epilink, "log" = exp, "identity" = identity)

    ## Clean the model environment when exiting the function
    on.exit(suppressWarnings(rm(cl, cumCIF, cumCIF.pb, data, doHessian,
        eventDists, eventsData, finetune, neghess, fisherinfo, fit, fixed,
        functions, globalEndemicIntercept, h.Intercept, inmfe, initpars,
        ll, negll, loglik, msgConvergence, msgNotConverged,
        mfe, mfhEvents, mfhGrid, model, my.na.action, na.action, namesOptimUser,
        namesOptimArgs, nlminbControl, nlminbRes, nlmObjective, nlmControl,
        nlmRes, nmRes, optim.args, optimArgs, control.siaf,
        optimMethod, optimRes, optimRes1, optimValid,
        origenv.endemic, origenv.epidemic, partial,
        partialloglik, ptm, qmatrix, res, negsc, score, start, subset, tmpexpr,
        typeSpecificEndemicIntercept, useScore, verbose, whichfixed, 
        inherits = FALSE)))
    
    ## also set fixed[st]iafpars to FALSE (for free posteriori evaluations, and
    ## to be defined for score function evaluation with optim.args=NULL)
    on.exit(fixedsiafpars <- fixedtiafpars <- FALSE, add = TRUE)


    ### Verify that 'data' inherits from "epidataCS"

    if (!inherits(data, "epidataCS")) {
        stop("'data' must inherit from class \"epidataCS\"")
    }


    ### Check time range

    if (!isScalar(t0) || !isScalar(T)) {
        stop("endpoints 't0' and 'T' must be single numbers")
    }
    if (T <= t0) {
        stop("'T' must be greater than 't0'")
    }
    if (!t0 %in% data$stgrid$start) {
        justBeforet0 <- match(TRUE, data$stgrid$start > t0) - 1L
        # if 't0' is beyond the time range covered by 'data$stgrid'
        if (is.na(justBeforet0)) justBeforet0 <- length(data$stgrid$start)   # t0 was too big
        if (justBeforet0 == 0L) justBeforet0 <- 1L   # t0 was too small
        t0 <- data$stgrid$start[justBeforet0]
        message("replaced 't0' by the value ", t0,
                " (must be a 'start' time of 'data$stgrid')")
    }
    if (!T %in% data$stgrid$stop) {
        justAfterT <- match(TRUE, data$stgrid$stop > T)
        # if 'T' is beyond the time range covered by 'data$stgrid'
        if (is.na(justAfterT)) justAfterT <- length(data$stgrid$stop)   # T was too big
        T <- data$stgrid$stop[justAfterT]
        message("replaced 'T' by the value ", T,
                " (must be a 'stop' time of 'data$stgrid')")
    }


    ### Subset events

    eventsData <- if (missing(subset)) data$events@data else {
        do.call("subset.data.frame", args = list(
            x = quote(data$events@data), subset = cl$subset, drop = FALSE
        ))
    }






    #############################################################
    ### Build up a model.frame for both components separately ###
    #############################################################



    ##########################
    ### epidemic component ###
    ##########################


    ### Parse epidemic formula

    if (missing(epidemic)) {
        origenv.epidemic <- parent.frame()
        epidemic <- ~ 0
    } else {
        origenv.epidemic <- environment(epidemic)
        environment(epidemic) <- environment()
        ## such that t0 and T are found in the subset expression below
    }
    epidemic <- terms(epidemic, data = eventsData, keep.order = TRUE)
    if (!is.null(attr(epidemic, "offset"))) {
        warning("offsets are not implemented for the 'epidemic' component")
    }


    ### Generate model frame

    # na.action mod such that for simulated epidataCS, where events of the
    # prehistory have missing 'BLOCK' indexes, those NA's do not matter.
    # ok because actually, 'eventBlocks' are only used in the partial likelihood
    # and there only eventBlocks[includes] is used (i.e. no prehistory events)
    my.na.action <- function (object, ...) {
        prehistevents <- row.names(object)[object[["(time)"]] <= t0]
        if (length(prehistevents) == 0L) return(na.action(object, ...))
        origprehistblocks <- object[prehistevents, "(BLOCK)"] # all NA
        object[prehistevents, "(BLOCK)"] <- 0L # temporary set non-NA
        xx <- na.action(object, ...)
        xx[match(prehistevents,row.names(xx),nomatch=0L), "(BLOCK)"] <-
            origprehistblocks[prehistevents %in% row.names(xx)]
        xx
    }

    mfe <- model.frame(epidemic, data = eventsData,
                       subset = time + eps.t > t0 & time <= T,
# here we can have some additional rows (individuals) compared to mfhEvents, which is established below!
# Namely those with time in (t0-eps.t; t0], i.e. still infective individuals, which are part of the prehistory of the process
                       na.action = my.na.action,
# since R 2.10.0 patched also works with epidemic = ~1 and na.action=na.fail (see PR#14066)
                       drop.unused.levels = FALSE,
                       time = time, tile = tile, type = type,
                       eps.t = eps.t, eps.s = eps.s, BLOCK = BLOCK,
                       obsInfLength = .obsInfLength, bdist = .bdist)

    
    ### Extract essential information from model frame

    # 'inmfe' indexes rows of data$events@data and is necessary for subsetting
    # influenceRegion (list incompatible with model.frame) and coordinates.
    # Note: model.frame() takes row.names from data
    inmfe <- which(row.names(data$events@data) %in% row.names(mfe))
    N <- length(inmfe)   # mfe also contains events of the prehistory
    eventTimes <- mfe[["(time)"]] # I don't use model.extract since it returns named vectors
    # Indicate events after t0, which are actually part of the process
    # (events in (-Inf;t0] only contribute in sum over infected individuals)
    includes <- which(eventTimes > t0)   # this indexes mfe!
    Nin <- length(includes)
    if (Nin == 0L) {
        stop("none of the ", nrow(data$events@data), " supplied ",
             "events is in the model (check 'subset', 't0' and 'T')")
    }
    eventBlocks <- mfe[["(BLOCK)"]]   # only necessary for partial log-likelihood
    eventTypes <- factor(mfe[["(type)"]])   # drop unused levels
    typeNames <- levels(eventTypes)
    nTypes <- length(typeNames)
    if (verbose && nTypes > 1L) cat("marked point pattern of", nTypes, "types\n")
    qmatrix <- checkQ(qmatrix, typeNames)
    # we only need the integer codes for the calculations
    eventTypes <- as.integer(eventTypes)


    ### Generate model matrix

    mme <- model.matrix(epidemic, mfe)
    q <- ncol(mme)
    hase <- q > 0L


    ### Extract further model components (only if q > 0)

    if (hase) {
        eps.t <- mfe[["(eps.t)"]]
        removalTimes <- eventTimes + eps.t
        eps.s <- mfe[["(eps.s)"]]
        bdist <- mfe[["(bdist)"]]
        gIntUpper <- mfe[["(obsInfLength)"]]
        gIntLower <- pmax(0, t0-eventTimes)
        eventCoords <- coordinates(data$events)[inmfe,,drop=FALSE]
        influenceRegion <- data$events@data$.influenceRegion[inmfe]
        iRareas <- vapply(X = influenceRegion, FUN = attr, which = "area",
                          FUN.VALUE = 0, USE.NAMES = FALSE)
        eventSources <- if (N == nobs(data) && identical(qmatrix, data$qmatrix)) {
            data$events@data$.sources
        } else { # re-determine because subsetting has invalidated row indexes
            ## if (verbose && N > 6000)
            ##     cat("calculating distance matrix of", N, "events ...\n")
            eventDists <- as.matrix(dist(eventCoords, method = "euclidean"))
            if (verbose) cat("updating list of potential sources ...\n")
            lapply(seq_len(N), function (i)
                determineSources(i, eventTimes, removalTimes, eventDists[i,],
                                 eps.s, eventTypes, qmatrix))
        }
        ## calculate sum_{k=1}^K q_{kappa_j,k} for all j = 1:N
        qSum <- unname(rowSums(qmatrix)[eventTypes])   # N-vector
    } else if (verbose) {
        message("no epidemic component in model")
    }


    ### Drop "terms" and restore original formula environment
    
    epidemic <- formula(epidemic)
    if (epilink != "log") # set as attribute only if non-standard link function
        attr(epidemic, "link") <- epilink
    environment(epidemic) <- origenv.epidemic
    ## We keep the original formula environment since it will be used to
    ## evaluate the modified twinstim-call in drop1/add1 (with default
    ## enclos=baseenv()), and cl$data should be visible from there.
    ## Alternatively, we could set it to parent.frame().




    #########################
    ### endemic component ###
    #########################


    ### Parse endemic formula

    if (missing(endemic)) {
        origenv.endemic <- parent.frame()
        endemic <- ~ 0
    } else {
        origenv.endemic <- environment(endemic)
        environment(endemic) <- environment()
        ## such that t0 and T are found in the subset expressions below
    }
    endemic <- terms(endemic, data = data$stgrid, keep.order = TRUE)

    ## check for type-specific endemic intercept and remove it from the formula
    ## (will be handled separately)
    typeSpecificEndemicIntercept <- "1 | type" %in% attr(endemic, "term.labels")
    if (typeSpecificEndemicIntercept) {
        endemic <- update(endemic, ~ . - (1|type)) # this drops the terms attributes
        endemic <- terms(endemic, data = data$stgrid, keep.order = TRUE)
    }

    globalEndemicIntercept <- if (typeSpecificEndemicIntercept) {
            attr(endemic, "intercept") <- 1L   # we need this to ensure that we have correct contrasts
            FALSE
        } else attr(endemic, "intercept") == 1L

    nbeta0 <- globalEndemicIntercept + typeSpecificEndemicIntercept * nTypes


    ### Generate endemic model frame and model matrix on event data

    mfhEvents <- model.frame(endemic, data = eventsData[row.names(mfe),],
                             subset = time>t0 & time<=T,
                             na.action = na.fail,
                             # since R 2.10.0 patched also works with
                             # endemic = ~1 (see PR#14066)
                             drop.unused.levels = FALSE)
    mmhEvents <- model.matrix(endemic, mfhEvents)
    # exclude intercept from endemic model matrix below, will be treated separately
    if (nbeta0 > 0) mmhEvents <- mmhEvents[,-1,drop=FALSE]
    #stopifnot(nrow(mmhEvents) == Nin)
    p <- ncol(mmhEvents)
    hash <- (nbeta0+p) > 0L


    ### Generate model frame and model matrix on grid data (only if p > 0)

    if (hash) {
        offsetEvents <- model.offset(mfhEvents)
        mfhGrid <- model.frame(endemic, data = data$stgrid,
                               subset = start >= t0 & stop <= T,
                               na.action = na.fail,
                               # since R 2.10.0 patched also works with
                               # endemic = ~1 (see PR#14066)
                               drop.unused.levels = FALSE,
                               BLOCK=BLOCK, tile=tile, dt=stop-start, ds=area)
                               # 'tile' is redundant here for fitting but useful
                               # for debugging & necessary for intensityplots
        gridBlocks <- mfhGrid[["(BLOCK)"]]
        histIntervals <- data$stgrid[!duplicated.default(
            data$stgrid$BLOCK, nmax = gridBlocks[length(gridBlocks)]
        ), c("BLOCK", "start", "stop")] # sorted
        row.names(histIntervals) <- NULL
        histIntervals <- histIntervals[histIntervals$start >= t0 &
                                       histIntervals$stop <= T,]
        gridTiles <- mfhGrid[["(tile)"]] # only needed for intensityplot
        mmhGrid <- model.matrix(endemic, mfhGrid)
        nGrid <- nrow(mmhGrid)

        # exclude intercept from endemic model matrix below, will be treated separately
        if (nbeta0 > 0) mmhGrid <- mmhGrid[,-1,drop=FALSE]
        # Extract endemic model components
        offsetGrid <- model.offset(mfhGrid)
        dt <- mfhGrid[["(dt)"]]
        ds <- mfhGrid[["(ds)"]]
        ## expression to calculate the endemic part on the grid -> .hIntTW()
        if (p > 0L) {
            hGridExpr <- quote(drop(mmhGrid %*% beta))
            if (!is.null(offsetGrid))
                hGridExpr <- call("+", quote(offsetGrid), hGridExpr)
        } else {
            hGridExpr <- if (is.null(offsetGrid))
                quote(numeric(nGrid)) else quote(offsetGrid)
        }
        hGridExpr <- call("exp", hGridExpr)
        ## expression to calculate the endemic part for the events -> .hEvents()
        hEventsExpr <- if (p > 0L) {
            quote(drop(mmhEvents %*% beta))
        } else {
            quote(numeric(Nin))
        }
        if (nbeta0 == 1L) { # global intercept
            hEventsExpr <- call("+", quote(beta0), hEventsExpr)
        } else if (nbeta0 > 1L) { # type-specific intercept
            hEventsExpr <- call("+", quote(beta0[eventTypes]), hEventsExpr)
        }
        if (!is.null(offsetEvents))
            hEventsExpr <- call("+", quote(offsetEvents), hEventsExpr)
        hEventsExpr <- call("exp", hEventsExpr)
    } else if (verbose) message("no endemic component in model")

    
    ### Drop "terms" and restore original formula environment
    
    endemic <- if (typeSpecificEndemicIntercept) {
        ## re-add it to the endemic formula
        update.formula(formula(endemic), ~ (1|type) + .)
    } else formula(endemic)
    environment(endemic) <- origenv.endemic
    ## We keep the original formula environment since it will be used to
    ## evaluate the modified twinstim-call in drop1/add1 (with default
    ## enclos=baseenv()), and cl$data should be visible from there.
    ## Alternatively, we could set it to parent.frame().


    ### Check that there is at least one parameter

    if (!hash && !hase) {
        stop("nothing to do: neither endemic nor epidemic parts were specified")
    }






    #############################
    ### Interaction functions ###
    #############################

    if (hase) {

        ## Check interaction functions
        siaf <- do.call(".parseiaf", args = alist(siaf, "siaf", eps.s, verbose))
        constantsiaf <- attr(siaf, "constant")
        nsiafpars <- siaf$npars

        tiaf <- do.call(".parseiaf", args = alist(tiaf, "tiaf", eps.t, verbose))
        constanttiaf <- attr(tiaf, "constant")
        ntiafpars <- tiaf$npars

        ## Check control.siaf
        if (constantsiaf) {
            control.siaf <- NULL
        } else if (is.list(control.siaf)) {
            if (!is.null(control.siaf$F)) stopifnot(is.list(control.siaf$F))
            if (!is.null(control.siaf$Deriv)) stopifnot(is.list(control.siaf$Deriv))
        } else if (!is.null(control.siaf)) {
            stop("'control.siaf' must be a list or NULL")
        }

        ## should we compute siafInt in parallel?
        useParallel <- cores > 1L && requireNamespace("parallel")
        ## but do not parallelize for a memoised siaf.step (becomes slower)
        if (useParallel &&
            !is.null(attr(siaf, "knots")) && !is.null(attr(siaf, "maxRange")) &&
            requireNamespace("memoise", quietly = TRUE) &&
            memoise::is.memoised(environment(siaf$f)$ringAreas)) {
            cores <- 1L
            useParallel <- FALSE
        }
        
        ## Define function that integrates the 'tiaf' function
        .tiafInt <- .tiafIntFUN()

        ## Define function that integrates the two-dimensional 'siaf' function
        ## over the influence regions of the events
        ..siafInt <- if (is.null(control.siaf[["siafInt"]])) {
            .siafInt <- .siafIntFUN(siaf = siaf, noCircularIR = all(eps.s > bdist),
                                    parallel = useParallel)
            ## Memoisation of .siafInt
            if (!constantsiaf && requireNamespace("memoise")) {
                memoise::memoise(.siafInt)
                ## => speed-up optimization since 'nlminb' evaluates the loglik and
                ## score for the same set of parameters at the end of each iteration
            } else {
                if (!constantsiaf && verbose)
                    message("Continuing without memoisation of 'siaf$f' cubature ...")
                .siafInt
            }
        } else {
            ## predefined cubature results in epitest(..., fixed = TRUE),
            ## where siafInt is identical during all permutations (only permuted)
            stopifnot(is.vector(control.siaf[["siafInt"]], mode = "numeric"),
                      length(control.siaf[["siafInt"]]) == N)
            local({
                env <- new.env(hash = FALSE, parent = .GlobalEnv)
                env$siafInt <- control.siaf[["siafInt"]]
                as.function(alist(siafpars=, ...=, siafInt), envir = env)
            })
        }
        .siafInt.args <- c(alist(siafpars), control.siaf$F)

    } else {
        
        if (!missing(siaf) && !is.null(siaf))
            warning("'siaf' can only be modelled in conjunction with an 'epidemic' process")
        if (!missing(tiaf) && !is.null(tiaf))
            warning("'tiaf' can only be modelled in conjunction with an 'epidemic' process")
        siaf <- tiaf <- NULL
        nsiafpars <- ntiafpars <- 0L
        control.siaf <- NULL

    }

    hassiafpars <- nsiafpars > 0L
    hastiafpars <- ntiafpars > 0L

    ## Can we calculate the score function?
    useScore <- if (partial) FALSE else if (hase) {
        (!hassiafpars | !is.null(siaf$deriv)) &
        (!hastiafpars | (!is.null(tiaf$deriv)) & !is.null(tiaf$Deriv))
    } else TRUE

    ## Define function that applies siaf$Deriv on all events (integrate the
    ## two-dimensional siaf$deriv function)
    if (useScore && hassiafpars) {
        .siafDeriv <- mapplyFUN(
            c(alist(siaf$Deriv, influenceRegion, type=eventTypes),
              list(MoreArgs=quote(list(siaf$deriv, siafpars, ...)),
                   SIMPLIFY=TRUE, USE.NAMES=FALSE)),
            ##<- we explicitly quote() the ...-part instead of simply including
            ##   it in the above alist() - only to make checkUsage() happy
            ## depending on nsiafpars, mapply() will return an N-vector
            ## or a nsiafpars x N matrix => transform to N x nsiafpars:
            after = quote(if (is.matrix(res)) t(res) else as.matrix(res)),
            parallel = useParallel)
        .siafDeriv.args <- c(alist(siafpars), control.siaf$Deriv)
    }



    ############################################################################
    ### Log-likelihood function, score function, expected Fisher information ###
    ############################################################################


    ### Total number of parameters (= length of 'theta')

    npars <- nbeta0 + p + q + nsiafpars + ntiafpars

    # REMINDER:
    #  theta - parameter vector c(beta0, beta, gamma, siafpars, tiafpars), where
    #    beta0   - endemic intercept (maybe type-specific)
    #    beta    - other parameters of the endemic component exp(offset + eta_h(t,s))
    #    gamma   - coefficients of the epidemic predictor
    #    siafpars- parameters of the epidemic spatial interaction function
    #    tiafpars- parameters of the epidemic temporal interaction function
    #  mmh[Events/Grid] - model matrix related to beta, i.e the endemic component,
    #                     either for events only or for the whole spatio-temporal grid
    #  offset[Events/Grid] - offset vector related to the endemic component (can be NULL),
    #                        either for events only or for the whole spatio-temporal grid
    #  dt, ds - columns of the spatio-temporal grid (dt = stop-start, ds = area)
    #  mme - model matrix related to gamma in the epidemic component
    #  siaf, tiaf - spatial/temporal interaction function (NULL, list or numeric)
    #  eventTimes, eventCoords, eventSources, gIntLower, gIntUpper, influenceRegion -
    #     columns of the events data frame


    if (hash)
    {
        ### Calculates the endemic component (for i in includes -> Nin-vector)
        ### h(t_i,s_i,kappa_i) = exp(offset_i + beta_{0,kappa_i} + eta_h(t_i,s_i))

        .hEvents <- function (beta0, beta) {}
        body(.hEvents) <- hEventsExpr


        ### Integral of the endemic component over [0;uppert] x W

        .hIntTW <- function (beta,
                             score = NULL, #matrix(1,nrow(mmhGrid),1L)
                             uppert = NULL) {}
        body(.hIntTW) <- as.call(c(as.name("{"),
            expression(
                subtimeidx <- if (!is.null(uppert)) { # && isScalar(uppert) && t0 <= uppert && uppert < T
                    if (uppert == t0) return(0)       # actually never happens
                                        # since uppert %in% eventTimes[includes] > t0
                    idx <- match(TRUE, histIntervals$stop >= uppert)
                    firstBlockBeyondUpper <- histIntervals$BLOCK[idx]
                    newdt <- uppert - histIntervals$start[idx]
                    dt[gridBlocks == firstBlockBeyondUpper] <- newdt
                    which(gridBlocks <= firstBlockBeyondUpper)
                } else NULL
            ),
            substitute(hGrid <- hGridExpr, list(hGridExpr=hGridExpr)),
            expression(sumterms <- hGrid * ds * dt),
            expression(if (is.null(score)) {
                if (is.null(subtimeidx))
                    sum(sumterms) else sum(sumterms[subtimeidx])
            } else {
                if (is.null(subtimeidx))
                    .colSums(score * sumterms, nGrid, ncol(score)) else
                .colSums((score * sumterms)[subtimeidx,,drop=FALSE], length(subtimeidx), ncol(score))
            })
        ))
    }

    if (hase)
    {
        ### Calculates the epidemic component for all events

        .eEvents <- function (gammapred, siafpars, tiafpars,
            ncolsRes = 1L, score = matrix(1,N,ncolsRes), f = siaf$f, g = tiaf$g)
            # second line arguments are for score functions with defaults for loglik
        {
            e <- vapply(X = includes, FUN = function (i) {
                sources <- eventSources[[i]]
                nsources <- length(sources)
                if (nsources == 0L) numeric(ncolsRes) else {
                    scoresources <- score[sources,,drop=FALSE]
                    predsources <- gammapred[sources]
                    repi <- rep.int(i, nsources)
                    sdiff <- eventCoords[repi,,drop=FALSE] - eventCoords[sources,,drop=FALSE]
                    fsources <- f(sdiff, siafpars, eventTypes[sources])
                    tdiff <- eventTimes[repi] - eventTimes[sources]
                    gsources <- g(tdiff, tiafpars, eventTypes[sources])
        # if(length(predsources) != NROW(fsources) || NROW(fsources) != NROW(gsources)) browser()
                    .colSums(scoresources * predsources * fsources * gsources,
                             nsources, ncolsRes)
                }
            }, FUN.VALUE = numeric(ncolsRes), USE.NAMES = FALSE)
            ## return a vector if ncolsRes=1, otherwise a matrix (Nin x ncolsRes)
            if (ncolsRes == 1L) e else t(e)
        }
    }


    ### Calculates the two components of the integrated intensity function
    ### over [0;uppert] x W x K

    heIntTWK <- function (beta0, beta, gammapred, siafpars, tiafpars,
                          uppert = NULL) {}
    body(heIntTWK) <- as.call(c(as.name("{"),
        if (hash) { # endemic component
            expression(
                hIntTW <- .hIntTW(beta, uppert = uppert),
                .beta0 <- rep_len(if (nbeta0==0L) 0 else beta0, nTypes),
                fact <- sum(exp(.beta0)),
                hInt <- fact * hIntTW
            )
        } else { expression(hInt <- 0) },
        if (hase) { # epidemic component
            c(expression(siafInt <- do.call("..siafInt", .siafInt.args)),#N-vector
              if (useParallel) expression( # print "try-catch"ed errors
                  if (any(.nonfinitesiafint <- !is.finite(siafInt)))
                  stop("invalid result of 'siaf$F' for 'siafpars=c(",
                       paste(signif(siafpars, getOption("digits")),
                             collapse=", "), ")':\n",
                       paste(unique(siafInt[.nonfinitesiafint]), sep="\n"),
                       call.=FALSE)
                  ),
              expression(
                  if (!is.null(uppert)) { # && isScalar(uppert) && t0 <= uppert && uppert < T
                      gIntUpper <- pmin(uppert-eventTimes, eps.t)
                      subtimeidx <- eventTimes < uppert
                      tiafIntSub <- .tiafInt(tiafpars,
                                             from = gIntLower[subtimeidx],
                                             to   = gIntUpper[subtimeidx],
                                             type = eventTypes[subtimeidx])
                      eInt <- sum(qSum[subtimeidx] * gammapred[subtimeidx] *
                                  siafInt[subtimeidx] * tiafIntSub)
                  } else {
                      tiafInt <- .tiafInt(tiafpars)
                      eInt <- sum(qSum * gammapred * siafInt * tiafInt)
                  }
                  )
              )
        } else expression(eInt <- 0),
        expression(c(hInt, eInt))
    ))


    ### Calculates the log-likelihood

    loglik <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        # dN part of the log-likelihood
        hEvents <- if (hash) .hEvents(beta0, beta) else 0
        eEvents <- if (hase) {
                gammapred <- drop(epilinkinv(mme %*% gamma)) # N-vector
                .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
            } else 0
        lambdaEvents <- hEvents + eEvents  # Nin-vector
        llEvents <- sum(log(lambdaEvents))
        # * llEvents is -Inf in case of 0-intensity at any event time
        # * If epilinkinv is 'identity', lambdaEvents < 0 if eEvents < -hEvents,
        #   and llEvents is NaN with a warning (intensity must be positive)
        if (is.nan(llEvents))  # nlminb() does not like NA function values
            llEvents <- -Inf

        # lambda integral of the log-likelihood
        heInt <- heIntTWK(beta0, beta, gammapred, siafpars, tiafpars)   # !hase => missing(gammapred), but lazy evaluation omits an error in this case because heIntTWK doesn't ask for gammapred
        llInt <- sum(heInt)

        # Return the log-likelihood
        ll <- llEvents - llInt
        ll
    }


    ### Calculates the score vector

    score <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        if (hase) {
            gammapred <- drop(epilinkinv(mme %*% gamma))  # N-vector
            hEvents <- if (hash) .hEvents(beta0, beta) else 0
            eEvents <- .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
            lambdaEvents <- hEvents + eEvents  # Nin-vector
            siafInt <- do.call("..siafInt", .siafInt.args) # N-vector
            tiafInt <- .tiafInt(tiafpars) # N-vector
        }

        # score vector for beta
        hScore <- if (hash)
        {
            score_beta0 <- if (nbeta0 == 1L) local({ # global intercept
                sEvents <- if (hase) {
                        hEvents / lambdaEvents
                    } else rep.int(1, Nin)
                sEventsSum <- sum(sEvents)
                sInt <- nTypes*exp(beta0) * .hIntTW(beta)
                sEventsSum - unname(sInt)
            }) else if (nbeta0 > 1L) local({ # type-specific intercepts
                ind <- sapply(seq_len(nTypes),
                              function (type) eventTypes == type,
                              simplify=TRUE, USE.NAMES=FALSE) # logical N x nTypes matrix
                sEvents <- if (hase) {
                        ind * hEvents / lambdaEvents
                    } else ind
                sEventsSum <- .colSums(sEvents, N, nTypes)
                sInt <- exp(beta0) * .hIntTW(beta)
                sEventsSum - unname(sInt)
            }) else numeric(0L) # i.e. nbeta0 == 0L

            score_beta <- if (p > 0L) local({
                sEvents <- if (hase) {
                        mmhEvents * hEvents / lambdaEvents
                    } else mmhEvents
                sEventsSum <- .colSums(sEvents, Nin, p)
                fact <- if (nbeta0 > 1L) sum(exp(beta0)) else if (nbeta0 == 1L) nTypes*exp(beta0) else nTypes
                sInt <- fact * .hIntTW(beta, mmhGrid)
                sEventsSum - sInt
            }) else numeric(0L)

            c(score_beta0, score_beta)
        } else numeric(0L)

        # score vector for gamma, siafpars and tiafpars
        eScore <- if (hase)
        {
            score_gamma <- local({
                nom <- .eEvents(switch(epilink, "log" = gammapred, "identity" = rep.int(1, N)),
                                siafpars, tiafpars,
                                ncolsRes=q, score=mme) # Nin-vector if q=1
                sEventsSum <- .colSums(nom / lambdaEvents, Nin, q)
                            # |-> dotted version also works for vector-arguments
                dgammapred <- switch(epilink, "log" = mme * gammapred, "identity" = mme)
                sInt <- .colSums(dgammapred * (qSum * siafInt * tiafInt), N, q)
                sEventsSum - sInt
            })

            score_siafpars <- if (hassiafpars && !fixedsiafpars) local({
                nom <- .eEvents(gammapred, siafpars, tiafpars,
                                ncolsRes=nsiafpars, f=siaf$deriv)
                sEventsSum <- .colSums(nom / lambdaEvents, Nin, nsiafpars)
                derivInt <- do.call(".siafDeriv", .siafDeriv.args) # N x nsiafpars matrix
                ## if useParallel, derivInt may contain "try-catch"ed errors
                ## in which case we receive a one-column character or list matrix
                if (!is.numeric(derivInt)) # we can throw a helpful error message
                    stop("invalid result of 'siaf$Deriv' for 'siafpars=c(",
                         paste(signif(siafpars, getOption("digits")),
                               collapse=", "), ")':\n",
                         paste(unique(derivInt[sapply(derivInt, is.character)]), sep="\n"),
                         call.=FALSE)
                sInt <- .colSums(derivInt * (qSum * gammapred * tiafInt),
                                 N, nsiafpars)
                sEventsSum - sInt
            }) else numeric(nsiafpars) # if 'fixedsiafpars', this part is unused

            score_tiafpars <- if (hastiafpars && !fixedtiafpars) local({
                nom <- .eEvents(gammapred, siafpars, tiafpars,
                                ncolsRes=ntiafpars, g=tiaf$deriv)
                sEventsSum <- .colSums(nom / lambdaEvents, Nin, ntiafpars)
                derivIntUpper <- tiaf$Deriv(gIntUpper, tiafpars, eventTypes)
                derivIntLower <- tiaf$Deriv(gIntLower, tiafpars, eventTypes)
                derivInt <- derivIntUpper - derivIntLower # N x ntiafpars matrix
                sInt <- .colSums(derivInt * (qSum * gammapred * siafInt), N, ntiafpars)
                sEventsSum - sInt
            }) else numeric(ntiafpars) # if 'fixedtiafpars', this part is unused

            c(score_gamma, score_siafpars, score_tiafpars)
        } else numeric(0L)

        # return the score vector
        scorevec <- c(hScore, eScore)
        scorevec
    }


    ### Estimates the expected Fisher information matrix
    ### by the "optional variation process" (Martinussen & Scheike, p. 64),
    ### or see Rathbun (1996, equation (4.7))

    fisherinfo <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        # only events (intdN) part of the score function needed
        zeromatrix <- matrix(0, Nin, 0)

        if (hase) {
            gammapred <- drop(epilinkinv(mme %*% gamma))  # N-vector
            hEvents <- if (hash) .hEvents(beta0, beta) else 0
            eEvents <- .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
            lambdaEvents <- hEvents + eEvents  # Nin-vector
        }

        # for beta
        hScoreEvents <- if (hash) {
            scoreEvents_beta0 <- if (nbeta0 > 1L) local({ # type-specific intercepts
                ind <- sapply(seq_len(nTypes),
                              function (type) eventTypes == type,
                              simplify=TRUE, USE.NAMES=FALSE) # logical N x nTypes matrix
                if (hase) {
                    ind * hEvents / lambdaEvents
                } else ind
            }) else if (nbeta0 == 1L) { # global intercept
                if (hase) {
                    hEvents / lambdaEvents
                } else matrix(1, Nin, 1L)
            } else zeromatrix

            scoreEvents_beta <- if (p > 0L) {
                if (hase) {
                    mmhEvents * hEvents / lambdaEvents
                } else mmhEvents   # Nin x p matrix
            } else zeromatrix

            unname(cbind(scoreEvents_beta0, scoreEvents_beta, deparse.level=0))
        } else zeromatrix

        # for gamma, siafpars and tiafpars
        eScoreEvents <- if (hase)
        {
            scoreEvents_gamma_nom <-
                .eEvents(switch(epilink, "log" = gammapred, "identity" = rep.int(1, N)),
                         siafpars, tiafpars, ncolsRes = q, score = mme)  # Ninxq matrix

            scoreEvents_siafpars_nom <- if (hassiafpars) {
                .eEvents(gammapred, siafpars, tiafpars, ncolsRes = nsiafpars, f = siaf$deriv)  # Ninxnsiafpars matrix
            } else zeromatrix

            scoreEvents_tiafpars_nom <- if (hastiafpars) {
                .eEvents(gammapred, siafpars, tiafpars, ncolsRes = ntiafpars, g = tiaf$deriv)  # Ninxntiafpars matrix
            } else zeromatrix

            eScoreEvents_nom <- cbind(scoreEvents_gamma_nom,
                                      scoreEvents_siafpars_nom,
                                      scoreEvents_tiafpars_nom,
                                      deparse.level=0)
            eScoreEvents_nom / lambdaEvents
        } else zeromatrix

        scoreEvents <- cbind(hScoreEvents, eScoreEvents, deparse.level=0)

        ## Build the optional variation process (Martinussen & Scheike, p64)
        ## info <- matrix(0, nrow = npars, ncol = npars,
        ##                dimnames = list(names(theta), names(theta)))
        ## for (i in 1:Nin) info <- info + crossprod(scoreEvents[i,,drop=FALSE])
        ## oh dear, this is nothing else but t(scoreEvents) %*% scoreEvents
        crossprod(scoreEvents)
    }


    ### Calculates the partial log-likelihood for continuous space
    ### (Diggle et al., 2009)

    partialloglik <- function (theta)
    {
        # Extract parameters from theta
        beta0    <- theta[seq_len(nbeta0)]
        beta     <- theta[nbeta0+seq_len(p)]
        gamma    <- theta[nbeta0+p+seq_len(q)]
        siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
        tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

        # calculcate the observed intensities
        hEvents <- if (hash) .hEvents(beta0, beta) else 0
        eEvents <- if (hase) {
                gammapred <- drop(epilinkinv(mme %*% gamma))  # N-vector
                .eEvents(gammapred, siafpars, tiafpars)  # Nin-vector! (only 'includes' here)
            } else 0
        lambdaEvents <- hEvents + eEvents  # Nin-vector

        # calculate integral of lambda(t_i, s, kappa) over at-risk set = (observation region x types)
        hInts <- if (hash) { # endemic component
                hGrid <- eval(hGridExpr)
                # integral over W and types for each time block in mfhGrid
                fact <- if (nbeta0 > 1L) sum(exp(beta0)) else if (nbeta0 == 1L) nTypes*exp(beta0) else nTypes
                hInt_blocks <- fact * tapply(hGrid*ds, gridBlocks, sum, simplify=TRUE)
                .idx <- match(eventBlocks[includes], names(hInt_blocks))
                unname(hInt_blocks[.idx])   # Nin-vector
            } else 0
        eInts <- if (hase) { # epidemic component
                siafInt <- do.call("..siafInt", .siafInt.args) # N-vector
                gs <- gammapred * siafInt # N-vector
                sapply(includes, function (i) {
                    timeSources <- determineSources(i, eventTimes, removalTimes,
                        0, Inf, NULL)
                    nSources <- length(timeSources)
                    if (nSources == 0L) 0 else {
                        repi <- rep.int(i, nSources)
                        tdiff <- eventTimes[repi] - eventTimes[timeSources]
                        gsources <- tiaf$g(tdiff, tiafpars, eventTypes[timeSources])
                        sum(qSum[timeSources] * gs[timeSources] * gsources)
                    }
                }, simplify=TRUE, USE.NAMES=FALSE)   # Nin-vector
            } else 0
        lambdaEventsIntW <- hInts + eInts   # Nin-vector

        # Calculate and return the partial log-likelihood
        p <- lambdaEvents / lambdaEventsIntW   # Nin-vector
        pll <- sum(log(p))
        pll
    }






    ################################
    ### Prepare for optimization ###
    ################################


    ll <- if (partial) partialloglik else loglik
    functions <- list(ll = ll,
                      sc = if (useScore) score else NULL,
                      fi = if (useScore) fisherinfo else NULL)

    
    ### Include check for validity of siafpars and tiafpars ('validpars') in ll

    if (!is.null(siaf$validpars)) {
        body(ll) <- as.call(append(as.list(body(ll)),
            as.list(expression(
                if (hassiafpars && !siaf$validpars(siafpars)) {
                    if (optimArgs$control$trace > 0L)
                        cat("(invalid 'siafpars' in loglik)\n")
                    return(-Inf)
                }
                )),
            after = grep("^siafpars <-", body(ll))))
    }

    if (!is.null(tiaf$validpars)) {
        body(ll) <- as.call(append(as.list(body(ll)),
            as.list(expression(
                if (hastiafpars && !tiaf$validpars(tiafpars)) {
                    if (optimArgs$control$trace > 0L)
                        cat("(invalid 'tiafpars' in loglik)\n")
                    return(-Inf)
                }
                )),
            after = grep("^tiafpars <-", body(ll))))
    }


    ### Check that optim.args is a list or NULL

    if (is.null(optim.args)) {  # no optimisation requested
        setting <- functions
        on.exit(rm(setting), add = TRUE)
        # Append model information
        setting$npars <- c(nbeta0 = nbeta0, p = p,
                           q = q, nsiafpars = nsiafpars, ntiafpars = ntiafpars)
        setting$qmatrix <- qmatrix   # -> information about nTypes and typeNames
        setting$formula <- list(endemic = endemic, epidemic = epidemic,
                                siaf = siaf, tiaf = tiaf)
        # Return settings
        setting$call <- cl
        environment(setting) <- environment()
        if (verbose)
            message("optimization skipped",
                    " (returning functions in data environment)")
        return(setting)
    } else if (!is.list(optim.args)) stop("'optim.args' must be a list or NULL")


    ### Check initial value for theta

    if (is.null(optim.args[["par"]])) { # set naive defaults
        h.Intercept <- if (nbeta0 > 0) rep.int(crudebeta0(
            nEvents = Nin,
            offset.mean = if (is.null(offsetGrid)) 0 else weighted.mean(offsetGrid, ds),
            W.area = sum(ds[gridBlocks==histIntervals[1,"BLOCK"]]),
            period = T-t0, nTypes = nTypes
            ), nbeta0) else numeric(0L)
        optim.args$par <- c(h.Intercept, rep.int(0, npars - nbeta0))
    } else { # check validity of par-specification
        if (!is.vector(optim.args$par, mode="numeric")) {
            stop("'optim.args$par' must be a numeric vector")
        }
        if (length(optim.args$par) != npars) {
            stop(gettextf(paste("'optim.args$par' (%d) does not have the same",
                                "length as the number of unknown parameters (%d)"),
                          length(optim.args$par), npars))
        }
    }
    
    ## Set names for theta
    names(optim.args$par) <- c(
        if (nbeta0 > 1L) {
            paste0("h.type",typeNames)
        } else if (nbeta0 == 1L) "h.(Intercept)",
        if (p > 0L) paste("h", colnames(mmhEvents), sep = "."),
        if (hase) paste("e", colnames(mme), sep = "."),
        if (hassiafpars) paste("e.siaf",1:nsiafpars,sep="."),
        if (hastiafpars) paste("e.tiaf",1:ntiafpars,sep=".")
    )

    ## values in "start" overwrite initial values given by optim.args$par
    if (!is.null(start)) {
        start <- check_twinstim_start(start)
        start <- start[names(start) %in% names(optim.args$par)]
        optim.args$par[names(start)] <- start
    }

    initpars <- optim.args$par


    ### Fixed parameters during optimization

    fixed <- optim.args[["fixed"]]
    optim.args[["fixed"]] <- NULL
    whichfixed <- if (is.null(fixed)) {
        integer(0L)
    } else if (isTRUE(fixed)) {
        seq_len(npars)
    } else {
        stopifnot(is.vector(fixed))
        if (is.numeric(fixed)) {
            stopifnot(fixed %in% seq_len(npars))
            fixed
        } else if (is.character(fixed)) {
            ## we silently ignore names of non-existent parameters
            intersect(fixed, names(initpars))
        } else if (is.logical(fixed)) {
            stopifnot(length(fixed) == npars)
            which(fixed)
        } else {
            stop("'optim.args$fixed' must be a numeric, character or logical vector")
        }
    }
    fixed <- setNames(logical(npars), names(initpars)) # FALSE
    fixed[whichfixed] <- TRUE
    fixedsiafpars <- hassiafpars && all(fixed[paste("e.siaf", 1:nsiafpars, sep=".")])
    fixedtiafpars <- hastiafpars && all(fixed[paste("e.tiaf", 1:ntiafpars, sep=".")])


    ### Define negative log-likelihood (score, hessian) for minimization
    ### as a function of the non-fixed parameters
    
    negll <- ll
    body(negll)[[length(body(negll))]] <-
        call("-", body(negll)[[length(body(negll))]])
    negsc <- if (useScore) {
        negsc <- score
        body(negsc)[[length(body(negsc))]] <-
            call("-", body(negsc)[[length(body(negsc))]])
        negsc
    } else NULL
    neghess <- if (useScore) fisherinfo else NULL

    if (any(fixed)) {
        ## modify negll, negsc and neghess for subvector optimization
        optim.args$par <- initpars[!fixed]
        if (verbose) {
            if (all(fixed)) {
                cat("\nno numerical likelihood optimization, all parameters fixed:\n")
            } else cat("\nfixed parameters during optimization:\n")
            print(initpars[fixed])
        }
        tmpexpr <- expression(
            initpars[!fixed] <- theta,
            theta <- initpars
            )
        body(negll) <- as.call(append(as.list(body(negll)), as.list(tmpexpr), 1))
        if (useScore) {
            body(negsc) <- as.call(append(as.list(body(negsc)), as.list(tmpexpr), 1))
            body(neghess) <- as.call(append(as.list(body(neghess)), as.list(tmpexpr), 1))
            # return non-fixed sub-vector / sub-matrix only
            body(negsc)[[length(body(negsc))]] <-
                call("[", body(negsc)[[length(body(negsc))]], quote(!fixed))
            body(neghess)[[length(body(neghess))]] <-
                call("[", body(neghess)[[length(body(neghess))]],
                     quote(!fixed), quote(!fixed), drop=FALSE)
        }

        ## if siafpars or tiafpars are fixed, pre-evaluate integrals    
        if (fixedsiafpars) {
            if (verbose)
                cat("pre-evaluating 'siaf' integrals with fixed parameters ...\n")
            if (!"memoise" %in% loadedNamespaces())
                cat("WARNING: Memoization of siaf integration not available!\n",
                    "         Repeated integrations with same parameters ",
                    "are redundant and slow!\n",
                    "         Really consider installing package \"memoise\"!\n",
                    sep="")
            siafInt <- local({
                siafpars <- initpars[paste("e.siaf", 1:nsiafpars, sep=".")]
                do.call("..siafInt", .siafInt.args) # memoise()d
            })
        }
        if (fixedtiafpars) {
            if (verbose) cat("pre-evaluating 'tiaf' integrals with fixed parameters ...\n")
            tiafInt <- .tiafInt(initpars[paste("e.tiaf", 1:ntiafpars, sep=".")])
            ## re-define .tiafInt such that it just returns the pre-evaluated
            ## integrals if called with the default arguments
            .tiafInt.orig <- .tiafInt
            body(.tiafInt) <- expression(
                if (nargs() == 1L) tiafInt else
                .tiafInt.orig(tiafpars, from, to, type, G)
            )
            ## restore the original function at the end
            on.exit({
                .tiafInt <- .tiafInt.orig
                rm(.tiafInt.orig)
            }, add=TRUE)
        }
    }






    if (any(!fixed)) {   ####################
                         ### Optimization ###
                         ####################

        ## Configure the optim procedure (check optim.args)

        # default arguments
        optimArgs <- list(par = NULL, # replaced by optim.args$par below
                          fn = quote(negll), gr = quote(negsc),
                          method = if (partial) "Nelder-Mead" else "nlminb",
                          lower = -Inf, upper = Inf,
                          control = list(), hessian = TRUE)
        # user arguments
        namesOptimArgs <- names(optimArgs)
        namesOptimUser <- names(optim.args)
        optimValid <- namesOptimUser %in% namesOptimArgs
        optimArgs[namesOptimUser[optimValid]] <- optim.args[optimValid]
        if (any(!optimValid)) {
            warning("unknown names in optim.args: ",
                    paste(namesOptimUser[!optimValid], collapse = ", "),
                    immediate. = TRUE)
        }
        doHessian <- optimArgs$hessian
        optimMethod <- optimArgs$method


        ## Call 'optim', 'nlminb', or 'nlm' with the above arguments

        if (verbose) {
            cat("\nminimizing the negative", if (partial) "partial", "log-likelihood",
                "using", if (optimMethod %in% c("nlm", "nlminb"))
                paste0("'",optimMethod,"()'") else {
                    paste0("'optim()'s \"", optimMethod, "\"")
                }, "...\n")
            cat("initial parameters:\n")
            print(optimArgs$par)
        }
        optimRes1 <- if (optimMethod == "nlminb") {
            nlminbControl <- control2nlminb(optimArgs$control,
                                            defaults = list(trace=1L, rel.tol=1e-6))
            ## sqrt(.Machine$double.eps) is the default reltol used in optim,
            ## which usually equals about 1.49e-08.
            ## The default rel.tol of nlminb (1e-10) seems too small
            ## (nlminb often does not finish despite no "relevant" change in loglik).
            ## I therefore use 1e-6, which is also the default in package nlme
            ## (see 'lmeControl').
            if (nlminbControl$trace > 0L) {
                cat("negative log-likelihood and parameters ")
                if (nlminbControl$trace == 1L) cat("in each iteration") else {
                    cat("every", nlminbControl$trace, "iterations") }
                cat(":\n")
            }
            nlminbRes <- nlminb(start = optimArgs$par, objective = negll,
                                gradient = negsc,
                                hessian = if (doHessian) neghess else NULL,
                                control = nlminbControl,
                                lower = optimArgs$lower, upper = optimArgs$upper)
            nlminbRes$value <- -nlminbRes$objective
            nlminbRes$counts <- nlminbRes$evaluations
            nlminbRes
        } else if (optimMethod == "nlm") {
            nlmObjective <- function (theta) {
                value <- negll(theta)
                grad <- negsc(theta)
                #hess <- neghess(theta)
                structure(value, gradient = grad)#, hessian = hess)
            }
            nlmControl <- optimArgs$control
            if (is.null(nlmControl[["print.level"]])) {
                nlmControl$print.level <- min(nlmControl$trace, 2L)
            }
            nlmControl$trace <- nlmControl$REPORT <- NULL
            if (is.null(nlmControl[["iterlim"]])) {
                nlmControl$iterlim <- nlmControl$maxit
            }
            nlmControl$maxit <- NULL
            nlmControl$check.analyticals <- FALSE
            ##<- we use the negative _expected_ Fisher information as the Hessian,
            ##   which is of course different from the true Hessian (=neg. obs. Fisher info)
            nlmRes <- do.call("nlm", c(alist(f = nlmObjective, p = optimArgs$par,
                                             hessian = doHessian),
                                             nlmControl))
            names(nlmRes)[names(nlmRes) == "estimate"] <- "par"
            nlmRes$value <- -nlmRes$minimum
            nlmRes$counts <- rep.int(nlmRes$iterations, 2L)
            nlmRes$convergence <- if (nlmRes$code %in% 1:2) 0L else nlmRes$code
            nlmRes
        } else { # use optim()
            optimArgs$control <- modifyList(list(trace=1L, REPORT=1L),
                                            optimArgs$control)
            if (finetune) optimArgs$hessian <- FALSE
            res <- do.call("optim", optimArgs)
            res$value <- -res$value
            res
        }


        ## Optional fine-tuning of ML estimates by robust Nelder-Mead

        optimRes <- if (finetune) {
            if (verbose) {
                cat("\nMLE from first optimization:\n")
                print(optimRes1$par)
                cat("loglik(MLE) =", optimRes1$value, "\n")
                cat("\nfine-tuning MLE using Nelder-Mead optimization ...\n")
            }
            optimArgs$par <- optimRes1$par
            optimArgs$method <- "Nelder-Mead"
            optimArgs$hessian <- doHessian
            optimArgs$control <- modifyList(list(trace=1L), optimArgs$control)
            nmRes <- do.call("optim", optimArgs)
            nmRes$value <- -nmRes$value
            nmRes$counts[2L] <- 0L   # 0 gradient evaluations (replace NA for addition below)
            nmRes
        } else optimRes1


        ## Convergence message
        
        msgConvergence <- if (finetune || optimMethod != "nlminb") {
            paste("code", optimRes$convergence)
        } else optimRes$message
        
        if (optimRes$convergence != 0) {
            msgNotConverged <- paste0("optimization routine did not converge (",
                                      msgConvergence, ")")
            warning(msgNotConverged)
            if (verbose) {
                cat("\nWARNING: ", msgNotConverged, "!\n", sep="")
                if ((finetune || optimMethod != "nlminb") &&
                    !is.null(optimRes$message) && nzchar(optimRes$message)) {
                    cat("MESSAGE: \"", optimRes$message, "\"\n", sep="")
                }
                if (hase && useScore && !constantsiaf &&
                    grepl("false", msgNotConverged)) {
                    cat("SOLUTION: increase the precision of 'siaf$Deriv' (and 'siaf$F')\n")
                    if (optimMethod == "nlminb") {
                        cat("          or nlminb's false convergence tolerance 'xf.tol'\n")
                    }
                }
            }
        }

        if (verbose) {
            cat("\n", if (finetune) "final ", "MLE:\n", sep = "")
            print(optimRes$par)
            cat("loglik(MLE) =", optimRes$value, "\n")
        }

    }






    ##############
    ### Return ###
    ##############


    ### Set up list object to be returned

    fit <- list(
           coefficients = if (any(fixed)) {
               if (all(fixed)) initpars else
               unlist(modifyList(as.list(initpars), as.list(optimRes$par)))
           } else optimRes$par,
           loglik = structure(if (all(fixed)) ll(initpars) else optimRes$value,
                              partial = partial),
           counts = if (all(fixed)) c("function"=1L, "gradient"=0L) else {
               optimRes1$counts + if (finetune) optimRes$counts else c(0L, 0L)
           },
           converged = if (all(fixed) || (optimRes$convergence == 0))
                       TRUE else msgConvergence
           )


    ### Add Fisher information matrices

    # estimation of the expected Fisher information matrix
    fit["fisherinfo"] <- list(
        if (useScore) structure(
            fisherinfo(fit$coefficients),
            dimnames = list(names(initpars), names(initpars))
            )
        )
    
    # If requested, add observed fisher info (= negative hessian at maximum)
    fit["fisherinfo.observed"] <- list(
        if (any(!fixed) && !is.null(optimRes$hessian)) optimRes$hessian
        ## no "-" here because we optimized the negative log-likelihood
        )


    ### Add fitted intensity values and integrated intensities at events

    # final coefficients
    theta    <- fit$coefficients
    beta0    <- theta[seq_len(nbeta0)]
    beta     <- theta[nbeta0+seq_len(p)]
    gamma    <- theta[nbeta0+p+seq_len(q)]
    siafpars <- theta[nbeta0+p+q+seq_len(nsiafpars)]
    tiafpars <- theta[nbeta0+p+q+nsiafpars+seq_len(ntiafpars)]

    # final siaf and tiaf integrals over influence regions / periods
    # and final gammapred (also used by intensity.twinstim)
    if (hase) {
        gammapred <- drop(epilinkinv(mme %*% gamma)) # N-vector
        if (!fixedsiafpars) siafInt <- do.call("..siafInt", .siafInt.args)
        if (!fixedtiafpars) tiafInt <- .tiafInt(tiafpars)
    }
    
    # fitted intensities
    hEvents <- if (hash) .hEvents(unname(beta0), beta) else rep.int(0, Nin)
    eEvents <- if (hase) {
            .eEvents(gammapred, siafpars, tiafpars) # Nin-vector! (only 'includes' here)
        } else rep.int(0, Nin)
    fit$fitted <- hEvents + eEvents   # = lambdaEvents  # Nin-vector
    fit$fittedComponents <- cbind(h = hEvents, e = eEvents)
    rm(hEvents, eEvents)
    
    # calculate cumulative ground intensities at event times
    # Note: this function is also used by residuals.twinstim
    LambdagEvents <- function (cores = 1L, cumCIF.pb = interactive())
    {
        if (cores != 1L) cumCIF.pb <- FALSE
        if (cumCIF.pb) pb <- txtProgressBar(min=0, max=Nin, initial=0, style=3)
        heIntEvents <- if (cores == 1L) {
            sapply(seq_len(Nin), function (i) {
                if (cumCIF.pb) setTxtProgressBar(pb, i)
                heIntTWK(beta0, beta, gammapred, siafpars, tiafpars,
                         eventTimes[includes[i]])
            }, simplify=TRUE, USE.NAMES=FALSE)
        } else {                        # cannot use progress bar
            simplify2array(parallel::mclapply(
                X=eventTimes[includes], FUN=heIntTWK,
                beta0=beta0, beta=beta,
                gammapred=gammapred, siafpars=siafpars,tiafpars=tiafpars,
                mc.preschedule=TRUE, mc.cores=cores
                ), higher=FALSE)
        }
        if (cumCIF.pb) close(pb)
        setNames(.colSums(heIntEvents, 2L, Nin), rownames(mmhEvents))
    }
    fit["tau"] <- list(
        if (cumCIF) {
            if (verbose)
                cat("\nCalculating fitted cumulative intensities at events ...\n")
            LambdagEvents(cores, cumCIF.pb)
        })

    # calculate observed R0's: mu_j = spatio-temporal integral of e_j(t,s) over
    # the observation domain (t0;T] x W (not whole R+ x R^2)
    fit$R0 <- if (hase) qSum * gammapred * siafInt * tiafInt else rep.int(0, N)
    names(fit$R0) <- row.names(mfe)

    
    ### Append model information

    fit$npars <- c(nbeta0 = nbeta0, p = p,
                   q = q, nsiafpars = nsiafpars, ntiafpars = ntiafpars)
    fit$qmatrix <- qmatrix   # -> information about nTypes and typeNames
    fit$bbox <- bbox(data$W)            # for completeness and for iafplot
    fit$timeRange <- c(t0, T)           # for simulate.twinstim's defaults
    fit$formula <- list(endemic = endemic, epidemic = epidemic,
                        siaf = siaf, tiaf = tiaf)
    fit["control.siaf"] <- list(control.siaf)    # might be NULL

    
    ### Append optimizer configuration

    optim.args$par <- initpars        # reset to also include fixed coefficients
    if (any(fixed)) optim.args$fixed <- names(initpars)[fixed] # restore
    fit$optim.args <- optim.args
    fit["functions"] <- list(
        if (model) {
            environment(fit) <- environment()
            functions
        })


    ### Return object of class "twinstim"

    if (verbose) cat("\nDone.\n")
    fit$call <- cl
    fit$runtime <- structure(proc.time() - ptm, cores=cores)
    class(fit) <- "twinstim"
    return(fit)

}
