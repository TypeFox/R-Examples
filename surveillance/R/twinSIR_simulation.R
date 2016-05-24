################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Simulate from a "twinSIR" model as described in Hoehle (2009)
###
### Copyright (C) 2009 Michael Hoehle, 2009, 2012, 2014 Sebastian Meyer
### $Revision: 1079 $
### $Date: 2014-10-18 01:26:00 +0200 (Sam, 18. Okt 2014) $
################################################################################

## Apart from simulation of SIR data, it is possible to simulate
## - SI: infPeriod = function(ids) rep(Inf, length(ids)
## - SIS: remPeriod = function(ids) rep(0, length(ids)
## - SIRS: remPeriod in (0;Inf)
##
## One can even simulate from a Cox model with the following settings:
## + no removal (i.e. infPeriod = function(ids) rep(Inf, length(ids))
## + no epidemic component (i.e. no alpha, no f, no w).

simEpidata <- function (formula, data, id.col, I0.col, coords.cols,
    subset, beta, h0, f = list(), w = list(), alpha, infPeriod,
    remPeriod = function(ids) rep(Inf, length(ids)),
    end = Inf, trace = FALSE, .allocate = NULL)
{
    cl <- match.call()
    
    #######################
    ### Check arguments ###
    #######################

    ### Build up model.frame
    mfnames <- c("", "formula", "data", "subset")
    mf <- cl[match(mfnames, names(cl), nomatch = 0L)]
    mf$na.action <- as.name("na.fail")
    mf$drop.unused.levels <- FALSE
    mf$xlev <- list()
    data <- eval(mf$data, parent.frame())
    if (!inherits(data, "data.frame")) {
        stop("'data' must inherit from class \"data.frame\"")
    }
    if (inherits(data, "epidata")) {
        id.col <- "id"
        I0.col <- "atRiskY"   # but we need !atRiskY (will be considered below)
        coords.cols <- names(data)[attr(data, "coords.cols")]
        if(length(formula) == 2L) { # i.e. no response specified
            formula[3L] <- formula[2L]
            formula[[2L]] <- quote(cbind(start, stop))
        }
    } else {
        for(colarg in c("id.col", "I0.col", "coords.cols")) {
            colidx <- get(colarg, inherits = FALSE)
            if (is.numeric(colidx)) {
                tmp <- names(data)[colidx]
                if (any(is.na(tmp))) {
                    stop("'", colarg, " = ", deparse(cl[[colarg]]), "': ", 
                         "column index must be in [1; ", ncol(data),
                         "=ncol(data)]")
                }
                assign(colarg, tmp, inherits = FALSE)
            }
        }
        mf$I0 <- if (is.null(I0.col)) {
                     substitute(rep(0, N), list(N=nrow(data)))
                 } else as.name(I0.col)
    }
    mf$id <- as.name(id.col)
    for(coords.col in coords.cols) {
        eval(call("$<-", quote(mf), coords.col, quote(as.name(coords.col))))
    }
    special <- c("cox")
    Terms <- terms(formula, specials = special, data = data,
                   keep.order = TRUE, simplify = FALSE)
    mf$formula <- Terms
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    ### Convert id to a factor (also removing unused levels if it was a factor)
    mf[["(id)"]] <- factor(mf[["(id)"]])
    ids <- levels(mf[["(id)"]])
    nObs <- length(ids)
    if (nObs == 0L) {
        stop("nothing to do: no individuals in 'data'")
    }
    idsInteger <- seq_len(nObs)
    
    ### Check start/stop consistency (response)
    .startstop <- model.response(mf)
    if (NCOL(.startstop) != 2L || !is.numeric(.startstop)) {
        stop("the lhs of 'formula' must be a numeric matrix with two columns ",
             "like 'cbind(start, stop)'")
    }
    timeIntervals <- unique(.startstop)
    timeIntervals <- timeIntervals[order(timeIntervals[,1L]), , drop = FALSE]
    nBlocks <- nrow(timeIntervals)
    if (any(timeIntervals[,2L] <= timeIntervals[,1L])) {
        stop("stop times must be greater than start times")
    }
    if (any(timeIntervals[-1L,1L] != timeIntervals[-nBlocks,2L])) {
        stop("inconsistent start/stop times: time intervals not consecutive")
    }
    
    ### Check .allocate
    if (is.null(.allocate)) {
        .allocate <- max(500, ceiling(nBlocks/100)*100)
    } else {
        if (!isScalar(.allocate) || .allocate < nBlocks) {
            stop("'.allocate' must be >= ", nBlocks)
        }
    }
    
    ### Check that all blocks are complete (all id's present)
    .blockidx <- match(.startstop[,1L], timeIntervals[,1L])
    if (any(table(.blockidx) != nObs)) {
        stop("all time intervals must be present for all id's")
    }
    
    ### Define a vector containing the time points where covariates change
    # unique 'start' time points (=> includes beginning of observation period)
    externalChangePoints <- as.vector(timeIntervals[,1L])
    
    ### SORT THE MODEL.FRAME BY BLOCK AND ID !!!
    mf <- mf[order(.blockidx, mf[["(id)"]]),]
    
    ### Extract the coordinates
    coords <- as.matrix(mf[idsInteger, tail(1:ncol(mf),length(coords.cols))])
    colnames(coords) <- coords.cols
    rownames(coords) <- ids
    
    ### Extract the endemic part Z of the design matrix (no intercept)
    des <- read.design(mf, Terms)
    Z <- des$Z
    nPredCox <- ncol(Z)   # number of endemic (cox) predictor terms
    
    ### Only include basic endemic variables in the event history output
    basicCoxNames <- rownames(attr(Terms,"factors"))[attr(Terms,"specials")$cox]
    basicVarNames <- sub("cox\\(([^)]+)\\)", "\\1", basicCoxNames)
    nBasicVars <- length(basicCoxNames)
    # this is necessary if some variables in 'formula' do not have main effects
    extraBasicVars <- as.matrix(mf[setdiff(basicCoxNames, colnames(Z))])
    
    ### Build up 3-dim array [id x time x var] of endemic terms 
    coxArray <- array(cbind(Z, extraBasicVars),
        dim = c(nObs, nBlocks, ncol(Z) + ncol(extraBasicVars)),
        dimnames = list(ids, NULL, c(colnames(Z), colnames(extraBasicVars))))
    idxPredVars <- seq_len(nPredCox)
    idxBasicVars <- match(basicCoxNames, dimnames(coxArray)[[3]])

    ### Check simulation parameters
    ## endemic (cox) part
    if (nPredCox > 0L) {
        if(missing(beta) || length(beta) != nPredCox || !is.numeric(beta)) {
            stop(gettextf(paste("'beta', a numeric vector of length %d",
                "(number of endemic terms), must be specified"), nPredCox))
        }
    } else {
        beta <- numeric(0L)
    }
    ## epidemic part
    nPredEpi <- length(f) + length(w)
    if (nPredEpi > 0L) {
        ## check f
        if (length(f) > 0L) {
            if (ncol(coords) == 0L) {
                stop("need coordinates for distance-based epidemic covariates 'f'")
            }
            if (!is.list(f) || is.null(names(f)) || any(!sapply(f, is.function))) {
                stop("'f' must be a named list of functions")
            }
            distmat <- as.matrix(dist(coords, method = "euclidean"))
        }
        ## check w
        if (length(w) > 0L) {
            if (!is.list(w) || is.null(names(w)) || any(!sapply(w, is.function))) {
                stop("'w' must be a named list of functions")
            }
            wijlist <- compute_wijlist(w = w, data = mf[idsInteger, ])
        }
        ## check alpha (coefficients for all of f and w)
        if (missing(alpha) || !is.numeric(alpha) || is.null(names(alpha))) {
            stop(gettextf(paste("'alpha', a named numeric vector of length %d",
                "(number of epidemic terms), must be specified"), nPredEpi))
        }
        alpha <- alpha[c(names(f), names(w))]
        if (any(is.na(alpha))) {
            stop("'alpha' is incomplete for 'f' or 'w'")
        }
        stopifnot(alpha >= 0)
    } else {
        alpha <- numeric(0L)
    }

    ### Parse the generator function for the infectious periods
    if (missing(infPeriod)) {
        stop("argument 'infPeriod' is missing (with no default)")
    }
    infPeriod <- match.fun(infPeriod)
    
    ### Parse the generator function for the removal periods
    remPeriod <- match.fun(remPeriod)

    ### Parse the log baseline function
    h0spec <- paste("'h0' must be a single number or a list of functions",
                    "\"exact\" and \"upper\"")
    if (missing(h0)) {
        stop(h0spec)
    }
    if (is.list(h0)) {
        if (!all(is.function(h0[["exact"]]), is.function(h0[["upper"]]))) {
            stop(h0spec)
        }
        if (!inherits(h0$upper, "stepfun")) {
            stop("function 'h0$upper' must be a 'stepfun'")
        }
        h0ChangePoints <- knots(h0$upper)
    } else if (isScalar(h0)) {
        h0func <- eval(parse(text = paste("function (t)", h0)))
        environment(h0func) <- parent.frame()
        h0 <- list(exact = h0func, upper = h0func)
        h0ChangePoints <- numeric(0L)
    } else {
        stop(h0spec)
    }
    if (!isScalar(h0$exact(0))) {
        stop("'h0$exact' must return a scalar")
    }
    
    ### Define function which decides if to reject a proposal during simulation
    exactEqualsUpper <- identical(h0$exact, h0$upper)
    mustReject <- if (exactEqualsUpper) {
            function () FALSE
        } else {
            function () lambdaStar/lambdaStarMax < runif(1)
        }

    ### Check simulation ending time
    if (!isScalar(end) || end <= 0) {
        stop("'end' must be a single positive numeric value")
    }


    ###################
    ### Preparation ###
    ###################
    
    ### Initialize set of infected and susceptible individuals
    infected <- which(
        mf[idsInteger,"(I0)"] == as.numeric(!inherits(data, "epidata"))
    ) # in case of "epidata", mf$(I0) equals data$atRiskY => infected = I0==0
    susceptibles <- which(! idsInteger %in% infected)
    
    ### Initialize tables of planned R-events and S-events
    Revents <- if (length(infected) > 0) {
        cbind(infected, infPeriod(ids[infected]))
    } else {
        matrix(numeric(0), ncol = 2)
    }
    Sevents <- matrix(numeric(0), ncol = 2)

    ### Small hook to subsequently update the (time depending) Cox predictor
    ### based on the current time (ct) during the simulation loop
    if (nPredCox > 0L) {
        coxUpdate <- expression(
            predCox <- as.matrix(
                    coxArray[,which(externalChangePoints == ct),idxPredVars]
                ) %*% beta
        )
    } else {
        predCox <- numeric(nObs)   # zeros
    }

    ### 'lambdaCalc' is the main expression for the calculation of the intensity 
    ### values IMMEDIATELY AFTER the current time 'ct'.
    ### It will be evaluated during the while-loop below.
    lambdaCalc <- expression(
        # Endemic Cox predictor (no h0 here!) of susceptibles
        predCoxS <- predCox[susceptibles],
        # Epidemic component of susceptibles
        lambdaEpidemic <- numeric(length(susceptibles)),   # zeros
        if (nPredEpi > 0L && length(infected) > 0L) {
            fCovars <- if (length(f) > 0L) {
                u <- distmat[,infected, drop = FALSE]
                vapply(X = f, FUN = function (B) rowSums(B(u)),
                       FUN.VALUE = numeric(nObs), USE.NAMES = FALSE)
            } else NULL
            wCovars <- if (length(w) > 0L) {
                vapply(X = wijlist, FUN = function (wij) {
                    rowSums(wij[, infected, drop = FALSE])
                }, FUN.VALUE = numeric(nobs), USE.NAMES = FALSE)
            } else NULL
            epiCovars <- cbind(fCovars, wCovars, deparse.level=0)
            # epiCovars is a matrix [nObs x nPredEpi] also used by updateNextEvent
            if (length(susceptibles) > 0L) {
                lambdaEpidemic <- epiCovars[susceptibles,,drop=FALSE] %*% alpha
            }
        },
        # Combined intensity
        lambdaS <- lambdaEpidemic + exp(h0$exact(ct) + predCoxS),
        # Ground intensity (sum of all lambdaS's)
        lambdaStar <- sum(lambdaS),
        # Upper bound on ground intensity
        lambdaStarMax <- if (exactEqualsUpper) {
                lambdaStar
            } else {
                sum(lambdaEpidemic) + sum(exp(h0$upper(ct) + predCoxS))
            }
    )
    # the following initializations are for R CMD check only ("visible binding")
    lambdaS <- numeric(length(susceptibles))
    lambdaStarMax <- lambdaStar <- numeric(1L)
    
    # At current time (ct) we have:
    # lambdaS is a _vector of length the current number of susceptibles_
    #     containing the intensity of infection for each susceptible individual.
    # lambdaStar is the overall infection rate.
    # lambdaStarMax is the upper bound for lambdaStar regarding baseline.
    # 'susceptible' and 'infected' are the corresponding sets of individuals
    #     immediately AFTER the last event
    # in theory, if a covariate changes in point t, then the intensity changes 
    # at t+0 only. intensities are left-continuous functions. time interval of
    # constant intensity is (start;stop]. but in the implementation we need at 
    # time ct the value of the log-baseline at ct+0, especially for 
    # ct %in% h0ChangePoints, thus h0$upper should be a stepfun with right=FALSE
    
    ### Create a history object alongside the simulation
    epiCovars0 <- matrix(0, nrow = nObs, ncol = nPredEpi,
                         dimnames = list(NULL, c(names(f), names(w))))
    basicVars0 <- matrix(0, nrow = nObs, ncol = nBasicVars,
                         dimnames = list(NULL, basicVarNames))
    emptyEvent <- cbind(BLOCK = 0, id = idsInteger,
                        start = 0, stop = 0, atRiskY = 0,
                        event = 0, Revent = 0, coords, basicVars0, epiCovars0)
    # WARNING: if you change the column order, you have to adjust the
    # hard coded column indexes everywhere below, also in getModel.simEpidata !
    .epiIdx <- tail(seq_len(ncol(emptyEvent)), nPredEpi)
    .basicIdx <- 7L + ncol(coords) + seq_len(nBasicVars)
    .nrowsEvHist <- .allocate * nObs   # initial size of the event history
    evHist <- matrix(NA_real_, nrow = .nrowsEvHist, ncol = ncol(emptyEvent),
                     dimnames = list(NULL, colnames(emptyEvent)))

    ## Hook - create new event and populate it with appropriate covariates
    updateNextEvent <- expression(
        nextEvent <- emptyEvent,
        # populate epidemic covariates
        if (nPredEpi > 0L && length(infected) > 0L) {
          nextEvent[,.epiIdx] <- epiCovars  # was calculated in lambdaCalc
        },
        # Which time is currently appropriate in (time varying) covariates
        tIdx <- match(TRUE, c(externalChangePoints,Inf) > ct) - 1L,
        if (nBasicVars > 0L) {
            nextEvent[,.basicIdx] <- coxArray[,tIdx,idxBasicVars]
        },
        # At-risk indicator
        if (length(susceptibles) > 0) {
            nextEvent[susceptibles,5L] <- 1
        },
        # Block index
        nextEvent[,1L] <- rep.int(block, nObs),
        # Start time
        nextEvent[,3L] <- rep.int(ct, nObs)
    )

    ## Hook function to add the event to the history
    addNextEvent <- expression(
        nextEvent[,4L] <- rep.int(ct, nObs), # stop time
        if (block*nObs > .nrowsEvHist) { # enlarge evHist if not big enough
            if (trace > 0L) {
                cat("Enlarging the event history @ block", block, "...\n")
            }
            evHist <- rbind(evHist,
              matrix(NA_real_, nrow = .allocate * nObs, ncol = ncol(emptyEvent))
            )
            .nrowsEvHist <- .nrowsEvHist + .allocate * nObs
        },
        evHistIdx <- idsInteger + nObs * (block-1),   # = seq.int(from = 1 + nObs*(block-1), to = nObs*block)
        evHist[evHistIdx,] <- nextEvent,
        block <- block + 1
    )


    #######################################################################
    ### MAIN PART: sequential simulation of infection and removal times ###
    #######################################################################
    
    ### Some indicators
    ct <- timeIntervals[1L,1L] # = externalChangePoints[1]   # current time
    block <- 1
    pointRejected <- FALSE
    loopCounter <- 0L
    trace <- as.integer(trace)
    hadNumericalProblemsInf <- hadNumericalProblems0 <- FALSE
    eventTimes <- numeric(0)

    ### Update (time depending) endemic covariates (if there are any)
    if (nPredCox > 0L) {
        eval(coxUpdate)
    }

    ### Let's rock 'n roll
    repeat {

        loopCounter <- loopCounter + 1L
        if (trace > 0L && loopCounter %% trace == 0L) {
            cat(loopCounter, "@t =", ct, ":\t|S| =", length(susceptibles),
                " |I| =", length(infected), "\trejected?", pointRejected, "\n")
        }

        if (!pointRejected) {
            ## Compute current conditional intensity
            eval(lambdaCalc)
            ## Update event history (uses epiCovars from lambdaCalc)
            eval(updateNextEvent)
        }
        pointRejected <- FALSE

        ## Determine time of next external change point
        changePoints <- c(externalChangePoints, h0ChangePoints,
                          Revents[,2], Sevents[,2])
        .isPendingChangePoint <- changePoints > ct
        nextChangePoint <- if (any(.isPendingChangePoint)) {
                min(changePoints[.isPendingChangePoint])
            } else Inf

        ## Simulate waiting time for the subsequent infection
        T <- tryCatch(rexp(1, rate = lambdaStarMax),
            warning = function(w) {
                if (!is.na(lambdaStarMax) && lambdaStarMax < 1) {   # rate was to small for rexp
                    if (length(susceptibles) > 0L) {
                        assign("hadNumericalProblems0", TRUE, inherits = TRUE)
                    }
                    if (nextChangePoint == Inf) NULL else Inf
                } else {   # rate was to big for rexp
                    0   # since R-2.7.0 rexp(1, Inf) returns 0 with no warning!
                }
            })

        ## Stop if lambdaStarMax too small AND no more changes in rate
        if (is.null(T)) {
            ct <- end
            eval(addNextEvent)
            break
        }
        
        ## Stop if lambdaStarMax too big meaning T == 0 (=> concurrent events)
        if (T == 0) {
            hadNumericalProblemsInf <- TRUE
            break
        }
        
        ## Stop at all costs if end of simulation time [0; end) has been reached
        if (isTRUE(min(ct+T, nextChangePoint) >= end)) {
        # ">=" because we don't want an event at "end"
            ct <- end
            eval(addNextEvent)
            break
        }

        if (ct + T > nextChangePoint) {
        ## Simulated time point is beyond the next time of intensity change
        ## (removal or covariate or upper baseline change point)
            ct <- nextChangePoint
            if (nPredCox > 0L && ct %in% externalChangePoints) {
                # update endemic covariates
                eval(coxUpdate)
            }
            if (.Reventidx <- match(ct, Revents[,2L], nomatch = 0L)) {
                # removal (I->R), thus update set of infected
                remover <- Revents[.Reventidx,1L]
                .remPeriod <- remPeriod(ids[remover])
                Sevents <- rbind(Sevents, c(remover, ct + .remPeriod))
                infected <- infected[-match(remover, infected)]
                nextEvent[remover,7L] <- 1
            }
            if (.Seventidx <- match(ct, Sevents[,2L], nomatch = 0L)) {
                # this will also be TRUE if above .remPeriod == 0 (SIS-like with pseudo-R-event)
                # re-susceptibility (R->S), thus update set of susceptibles
                resusceptible <- Sevents[.Seventidx,1L]
                susceptibles <- c(susceptibles, resusceptible)
            }
            # update event history
            eval(addNextEvent)
        } else {
        ## Simulated time point lies within the thinning period
        ## => rejection sampling step
            ct <- ct + T
            if (length(h0ChangePoints) > 0L) {# i.e. if non-constant baseline
                # Update intensities for rejection probability at new ct
                eval(lambdaCalc)
            }
            if (mustReject()) {
                pointRejected <- TRUE
                next
            }
            # At this point, we have an actual event! =>
            # Sample the individual who becomes infected with probabilities
            # according to the intensity proportions
            victimSindex <- sample(length(susceptibles), 1L,
                                   prob = lambdaS/lambdaStar)
            victim <- susceptibles[victimSindex]
            eventTimes <- c(eventTimes, ct)
            Revents <- rbind(Revents, c(victim, ct + infPeriod(ids[victim])))
            susceptibles <- susceptibles[-victimSindex]
            infected <- c(infected, victim)
            
            # Add to history
            nextEvent[victim,6L] <- 1
            eval(addNextEvent)
        }
    
    }


    ##############
    ### Return ###
    ##############
    
    if (hadNumericalProblemsInf) {
        warning("simulation ended due to an infinite overall infection rate")
    }
    if (hadNumericalProblems0) {
        warning("occasionally, the overall infection rate was numerically ",
                "equal to 0 although there were individuals at risk")
    }
    
    if (trace > 0L) {
      cat("Converting the event history into a data.frame (\"epidata\") ...\n")
    }
    epi <- as.data.frame(evHist[seq_len(nObs*(block-1)),,drop=FALSE])
    epi$id <- factor(ids[epi$id], levels = ids)
    rownames(epi) <- NULL
    attr(epi, "eventTimes") <- eventTimes
    attr(epi, "timeRange") <- c(timeIntervals[1L,1L], ct)
    attr(epi, "coords.cols") <- 7L + seq_len(ncol(coords))
    attr(epi, "f") <- f
    attr(epi, "w") <- w
    attr(epi, "config") <- list(h0 = h0$exact, beta = beta, alpha = alpha)
    attr(epi, "call") <- cl
    attr(epi, "terms") <- Terms
    class(epi) <- c("simEpidata", "epidata", "data.frame")
    if (trace > 0L) {
        cat("Done.\n")
    }
    return(epi)
}


### We define no plot-method for simEpidata (as a wrapper for intensityPlot),
### because we want plot(simEpidataObject) to use the inherited method plot.epidata
### which shows the evolution of the numbers of individuals in states S, I, and R


################################################################################
# A 'simulate' method for objects of class "twinSIR".
################################################################################

simulate.twinSIR <- function (object, nsim = 1, seed = 1,
    infPeriod = NULL, remPeriod = NULL,
    end = diff(range(object$intervals)), trace = FALSE, .allocate = NULL,
    data = object$data, ...)
{
    theta <- coef(object)
    px <- ncol(object$model$X)
    pz <- ncol(object$model$Z)
    nh0 <- attr(object$terms, "intercept") * length(object$nEvents)
    
    f <- object$model$f  # contains only the f's used in the model formula
    w <- object$model$w  # contains only the w's used in the model formula
    if (any(missingf <- !names(f) %in% colnames(object$model$X))) {
        stop("simulation requires distance functions 'f', missing for: ",
             paste(colnames(object$model$X)[missingf], collapse=", "))
    }
    if (any(missingw <- !names(w) %in% colnames(object$model$X))) {
        stop("simulation requires functions 'w', missing for: ",
             paste(colnames(object$model$X)[missingw], collapse=", "))
    }

    formulaLHS <- "cbind(start, stop)"
    formulaRHS <- paste(c(as.integer(nh0 > 0), # endemic intercept?
                          names(theta)[px+nh0+seq_len(pz-nh0)]),
                        collapse = " + ")
    formula <- formula(paste(formulaLHS, formulaRHS, sep="~"),
                       env = environment(formula(object)))
    h0 <- if (nh0 == 0L) {
              if (pz == 0L) {
                -Inf   # no endemic component at all (exp(-Inf) == 0)
              } else {
                0      # endemic covariates act on 0-log-baseline hazard
              }
          } else {
              .h0 <- stepfun(x = object$intervals[1:nh0],
                             y = c(0,theta[px+seq_len(nh0)]),
                             right = FALSE)
              list(exact = .h0, upper = .h0)
          }
    
    if (!inherits(data, "epidata")) {
        stop("invalid 'data' argument: use function 'twinSIR' with ",
             "'keep.data = TRUE'")
    }
    if (is.null(infPeriod) || is.null(remPeriod)) {
        s <- summary(data)
        eventsByID <- s$byID
        if (is.null(infPeriod)) {
            infPeriod <-
            if (s$type == "SI") {
                function (ids) rep.int(Inf, length(ids))
            } else { # SIR, SIRS or SIS
                eventsByID$infPeriod <- eventsByID$time.R - eventsByID$time.I
                meanInfPeriodByID <- if (s$type %in% c("SIRS", "SIS")) {
                        c(tapply(eventsByID$infPeriod, list(eventsByID$id),
                                 mean, na.rm = TRUE, simplify = TRUE))
                    } else {
                        structure(eventsByID$infPeriod, names = eventsByID$id)
                    }
                meanInfPeriod <- mean(meanInfPeriodByID, na.rm = TRUE)
                if (is.na(meanInfPeriod)) {
                    stop("'infPeriod = NULL' invalid: ",
                         "no infection periods observed")
                }
                function (ids) {
                    infPeriods <- meanInfPeriodByID   # named vector
                    infPeriods[is.na(infPeriods)] <- meanInfPeriod
                    infPeriods[ids]
                }
            }
        }
        if (is.null(remPeriod)) {
            remPeriod <-
            if (s$type == "SIRS") {
                eventsByID$remPeriod <- eventsByID$time.S - eventsByID$time.R
                meanRemPeriodByID <- c(tapply(eventsByID$remPeriod,
                    list(eventsByID$id), mean, na.rm = TRUE, simplify = TRUE))
                meanRemPeriod <- mean(meanRemPeriodByID, na.rm = TRUE)
                function (ids) {
                    remPeriods <- meanRemPeriodByID   # named vector
                    remPeriods[is.na(remPeriods)] <- meanRemPeriod
                    remPeriods[ids]
                }
            } else if (s$type == "SIS") {
                function (ids) rep.int(0, length(ids))
            } else { # SIR or SI
                function (ids) rep.int(Inf, length(ids))
            }
        }
    }
    set.seed(seed)
    res <- replicate(nsim,
        simEpidata(formula, data = data,
                   beta = theta[px + nh0 + seq_len(pz-nh0)],
                   h0 = h0, f = f, w = w, alpha = theta[seq_len(px)],
                   infPeriod = infPeriod, remPeriod = remPeriod,
                   end = end, trace = trace, .allocate = .allocate),
        simplify = FALSE
    )
    if (nsim == 1L) res[[1L]] else res
}

