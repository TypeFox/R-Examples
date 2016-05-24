################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Simulate a point pattern according to a spatio-temporal intensity model of
### class "twinstim". The function basically uses Ogata's modified thinning
### algorithm (cf. Daley & Vere-Jones, 2003, Algorithm 7.5.V.).
###
### Copyright (C) 2010-2016 Sebastian Meyer
### $Revision: 1608 $
### $Date: 2016-03-04 22:40:13 +0100 (Fre, 04. Mär 2016) $
################################################################################

### CAVE:
### - the type of contrasts for factor variables has to be set through options("contrasts")
### - if epidemic-only process (!hash), we actually don't need stgrid, but we
###   want to have valid epidataCS at the end, which requires stgrid

## model.frame() evaluates '...' with 'data'
utils::globalVariables(c("BLOCK", "tile", "area"))

simEpidataCS <- function (endemic, epidemic, siaf, tiaf, qmatrix, rmarks,
    events, stgrid, tiles, beta0, beta, gamma, siafpars, tiafpars, epilink = "log",
    t0 = stgrid$start[1], T = tail(stgrid$stop,1), nEvents = 1e5,
    control.siaf = list(F=list(), Deriv=list()),
    W = NULL, trace = 5, nCircle2Poly = 32, gmax = NULL, .allocate = 500,
    .skipChecks = FALSE, .onlyEvents = FALSE)
{
    ptm <- proc.time()[[3]]
    cl <- match.call()


    #######################
    ### Check arguments ### (this takes many lines of code ...)
    #######################


    cat("\nChecking the supplied arguments ...\n")

    ### Some simple input checks

    if (missing(endemic)) endemic <- ~ 0 else stopifnot(inherits(endemic, "formula"))
    if (missing(epidemic)) epidemic <- ~ 0 else stopifnot(inherits(epidemic, "formula"))
    if (length(trace) != 1L) stop("'trace' must be a single integer or logical value")
    trace <- as.integer(trace)
    if (!isScalar(nCircle2Poly)) stop("'nCircle2Poly' must be scalar")
    nCircle2Poly <- as.integer(nCircle2Poly)
    if (!isScalar(.allocate)) stop("'.allocate' must be scalar")
    .allocate <- as.integer(.allocate)
    .skipChecks <- as.logical(.skipChecks)
    .onlyEvents <- as.logical(.onlyEvents)

    
    ### Check qmatrix

    if (missing(qmatrix)) qmatrix <- diag(1)
    nTypes <- nrow(qmatrix)
    if (is.null(typeNames <- rownames(qmatrix))) {
        if (nTypes > length(LETTERS)) stop("'qmatrix' needs dimnames")
        typeNames <- LETTERS[seq_len(nTypes)]
    }
    qmatrix <- checkQ(qmatrix, typeNames)
    qSumTypes <- rowSums(qmatrix)  # how many types can be triggered by each type


    ### Check other "epidataCS" components (events, stgrid, tiles, and W)
    
    if (!missing(events) && !is.null(events)) {
        events <- events[!names(events) %in% reservedColsNames_events]
        if (!.skipChecks) {
            cat("Checking 'events':\n")
            events <- check_events(events, dropTypes = FALSE)
            # epscols are obligatory in 'check_events', which is also appropriate here
        }
        ## check event types
        events@data$type <- factor(events@data$type, levels=typeNames)
        if (any(.typeIsNA <- is.na(events@data$type))) {
            warning("ignored some 'events' of unknown type")
            events <- events[!.typeIsNA,]
        }
    }

    if (!.skipChecks) {
        cat("Checking 'stgrid':\n")
        stgrid <- check_stgrid(stgrid[grep("^BLOCK$", names(stgrid), invert=TRUE)])
    }

    W <- if (is.null(W)) {
        cat("Building 'W' as the union of 'tiles' ...\n")
    	unionSpatialPolygons(tiles)
    } else check_W(W)  # does as(W, "SpatialPolygons")

    tileLevels <- levels(stgrid$tile)
    tiles <- check_tiles(tiles, tileLevels,
                         areas.stgrid = stgrid[["area"]][seq_along(tileLevels)],
                         W = W, keep.data = FALSE)
    
    ## Transform W to class "owin"
    Wowin <- as(W, "owin")
    Wedges <- edges(Wowin, check = FALSE)
    maxExtentOfW <- diameter.owin(Wowin)

    
    ### Check parameters

    beta0 <- if (missing(beta0)) numeric(0L) else as.vector(beta0, mode="numeric")
    beta  <- if (missing(beta))  numeric(0L) else as.vector(beta,  mode="numeric")
    gamma <- if (missing(gamma)) numeric(0L) else as.vector(gamma, mode="numeric")
    siafpars <- if (missing(siafpars)) numeric(0L) else as.vector(siafpars, mode="numeric")
    tiafpars <- if (missing(tiafpars)) numeric(0L) else as.vector(tiafpars, mode="numeric")
    nbeta0 <- length(beta0)
    if (nbeta0 > 1L && nbeta0 != nTypes) {
        stop("'beta0' must have length 0, 1, or 'nrow(qmatrix)'")
    }
    p <- length(beta)
    q <- length(gamma)
    nsiafpars <- length(siafpars)
    ntiafpars <- length(tiafpars)

    hase <- q > 0L
    hassiafpars <- nsiafpars > 0L
    hastiafpars <- ntiafpars > 0L
    if (!hase && (hassiafpars | hastiafpars)) {
        stop("'siafpars' and 'tiafpars' require 'gamma'")
    }


    ### Check time range

    if (is.null(t0)) t0 <- eval(formals()$t0)
    if (is.null(T)) T <- eval(formals()$T)
    if (!isScalar(t0) || !isScalar(T)) {
        stop("endpoints 't0' and 'T' must be single numbers")
    }
    if (T <= t0) {
        stop("'T' must be greater than 't0'")
    }
    stopifnot(t0 >= stgrid$start[1], T <= tail(stgrid$stop,1))


    ### Subset stgrid to include actual time range only

    # BLOCK in stgrid such that start time is equal to or just before t0
    block_t0 <- stgrid$BLOCK[match(TRUE, c(stgrid$start,Inf) > t0) - 1L]
    # BLOCK in stgrid such that stop time is equal to or just after T
    block_T <- stgrid$BLOCK[match(TRUE, stgrid$stop >= T)]
    stgrid <- stgrid[stgrid$BLOCK>=block_t0 & stgrid$BLOCK<=block_T,,drop=FALSE]
    stgrid$start[stgrid$BLOCK == block_t0] <- t0
    stgrid$stop[stgrid$BLOCK == block_T] <- T
    # matrix of BLOCKS and start times (used later)
    blockstarts <- with(stgrid,
        cbind(block_t0:block_T, start[match(block_t0:block_T, BLOCK)],
              deparse.level = 0L)
    )


    ### Check mark-generating function

    # eps.t and eps.s are also unpredictable marks (generated by rmarks)
    unpredMarks <- unique(c("eps.t", "eps.s", if (hase) {
        setdiff(all.vars(epidemic), c("type", names(stgrid)))
    }))
    rmarks <- match.fun(rmarks)
    sampleCoordinate <- coordinates(spsample(tiles, n=1L, type="random"))
    sampleMarks <- rmarks(t0, sampleCoordinate)
    # should be a one-row data.frame
    if (!is.data.frame(sampleMarks) || nrow(sampleMarks) != 1L) {
        stop("'rmarks' must return a one-row data.frame of marks")
    }
    markNames <- names(sampleMarks)
    if (.idx <- match(FALSE, unpredMarks %in% markNames, nomatch=0L)) {
        stop("the unpredictable mark '", unpredMarks[.idx], "' is not returned by 'rmarks'")
    }
    if (!all(sapply(sampleMarks[unpredMarks], function(x)
                    inherits(x, c("integer","numeric","logical","factor"), which=FALSE))))
        warning("'rmarks' should return \"numeric\", \"logical\", or",
                " \"factor\" ('epidemic') variables only")


    ### Check prehistory of the process

    Nout <- 0L
    if (!missing(events) && !is.null(events)) {
        .stillInfective <- with(events@data, time <= t0 & time + eps.t > t0)
        Nout <- sum(.stillInfective)
        events <- if (Nout > 0L) {
            events[.stillInfective,]
        } else {
            .eventstxt <- if (.skipChecks) "data$events" else "events"  # for simulate.twinstim
            cat("(no events from '", .eventstxt,
                "' were considered as prehistory)\n", sep="")
            NULL
        }
    }
    
    ## separate coordinates and data
    if (Nout > 0L) {
        check_tiles_events(tiles, events)
        eventCoords <- coordinates(events)
        eventData <- events@data
        ## check presence of unpredictable marks
        if (length(.idx <- which(!unpredMarks %in% names(eventData)))) {
            stop("missing unpredictable marks in 'events': ",
                 paste0("\"", unpredMarks[.idx], "\"", collapse=", "))
        }
        ## check type of unpredictable marks
        for (um in unpredMarks) {
            if (!identical(class(sampleMarks[[um]]), class(eventData[[um]])))
                stop("the class of the unpredictable mark '", um,
                     "' in the 'events' prehistory ",
                     "is not identical to the class returned by 'rmarks'")
        }
        ## add marks which are not in the prehistory but simulated by 'rmarks'
        if (length(.add2events <- setdiff(markNames, names(eventData)))) {
            eventData <- cbind(eventData, sampleMarks[.add2events])
            is.na(eventData[.add2events]) <- TRUE
        }
        eventData <- eventData[c("time", "tile", "type", markNames)]
    } else { ## empty prehistory
        eventCoords <- matrix(0, nrow=0L, ncol=2L)
        eventData <- data.frame(
            time = numeric(0L),
            tile = factor(character(0L), levels=tileLevels),
            type = factor(character(0L), levels=typeNames),
            check.rows = FALSE, check.names = FALSE
            )
        eventData <- cbind(eventData, sampleMarks[0L,])
    }

    ## helper function to attach covariates from 'stgrid' to events
    attachstgridvars <- function (eventData, stgridvars)
    {
        if (length(stgridvars) == 0L) return(eventData)
        gridcellsOfEvents <- integer(nrow(eventData))
        for (i in seq_along(gridcellsOfEvents)) {
            gridcellsOfEvents[i] <- gridcellOfEvent(eventData[i,"time"],
                                                    eventData[i,"tile"],
                                                    stgrid)
        }
        cbind(eventData, stgrid[gridcellsOfEvents, stgridvars, drop=FALSE])
    }


    ### Build epidemic model matrix

    epidemic <- terms(epidemic, data = eventData, keep.order = TRUE)
    if (!is.null(attr(epidemic, "offset"))) {
        warning("offsets are not implemented for the 'epidemic' component")
    }

    # helper function taking eventData and returning the epidemic model.matrix
    buildmme <- function (eventData)
    {
        # which variables do we have to copy from stgrid?
        stgridCopyCols <- match(all.vars(epidemic), names(stgrid), nomatch = 0L)
        eventData <- attachstgridvars(eventData, stgridCopyCols)
        mfe <- model.frame(epidemic, data = eventData,
                           na.action = na.fail, drop.unused.levels = FALSE)
        model.matrix(epidemic, mfe)
    }
    mme <- buildmme(eventData)

    if (ncol(mme) != q) {
        cat(ncol(mme), "epidemic model terms:\t",
            paste(colnames(mme), collapse="  "), "\n")
        stop("length of 'gamma' (", q,
             ") does not match the 'epidemic' specification (",
             ncol(mme), ")")
    }

    ## (inverse) link function for the epidemic linear predictor of event marks
    epilink <- match.arg(epilink, choices = c("log", "identity"))
    epilinkinv <- switch(epilink, "log" = exp, "identity" = identity)


    ### Build endemic model matrix

    endemic <- terms(endemic, data = stgrid, keep.order = TRUE)

    # check if we have an endemic component at all
    hasOffset <- !is.null(attr(endemic, "offset"))
    hash <- (nbeta0 + p + hasOffset) > 0L
    if (!hash) {
        if (!hase) {
            stop("nothing to do: neither endemic nor epidemic parameters were specified")
            # actually, the process might be endemic offset-only, which I don't care about ATM
        }
        if (Nout == 0L) {
            stop("missing 'events' pre-history (no endemic component)")
        }
    }

    # remove (1|type) specification
    typeSpecificEndemicIntercept <-
        "1 | type" %in% attr(endemic, "term.labels") || nbeta0 > 1
    if (typeSpecificEndemicIntercept) {
        endemic <- update(endemic, ~ . - (1|type)) # this drops the terms attributes
        endemic <- terms(endemic, data = stgrid, keep.order = TRUE)
        if (nbeta0 <= 1L) {
            stop("for type-specific endemic intercepts, 'beta0' must be longer than 1")
        }
    }

    # ensure that we have correct contrasts in the endemic component
    attr(endemic, "intercept") <- as.integer(nbeta0 > 0L)

    # helper function taking eventData (with time and tile columns)
    # and returning the endemic model.matrix
    buildmmh <- function (eventData)
    {
        # if 'pi' appears in 'endemic' we don't care, and if a true covariate is
        # missing, model.frame will throw an error
        # which variables do we have to copy from stgrid?
        stgridCopyCols <- match(all.vars(endemic), names(stgrid), nomatch = 0L)

        # attaching covariates from 'stgrid' to events
        eventData <- attachstgridvars(eventData, stgridCopyCols)
        # construct model matrix
        mfhEvents <- model.frame(endemic, data = eventData, na.action = na.fail,
                                 drop.unused.levels = FALSE)
        mmhEvents <- model.matrix(endemic, mfhEvents)
        # exclude intercept from endemic model matrix below, will be treated separately
        if (nbeta0 > 0) mmhEvents <- mmhEvents[,-1,drop=FALSE]
        structure(mmhEvents, offset = model.offset(mfhEvents))
    }

    # actually, we don't need the endemic model matrix for the pre-history events at all
    # this is just to test consistence with 'beta' and for the names of 'beta'
    mmh <- buildmmh(eventData[0L,])
    if (ncol(mmh) != p) {
        stop("length of 'beta' (", p, ") does not match the 'endemic' specification (", ncol(mmh), ")")
    }


    ### Build endemic model matrix on stgrid

    mfhGrid <- model.frame(endemic, data = stgrid, na.action = na.fail,
                           drop.unused.levels = FALSE,
                           BLOCK = BLOCK, tile = tile, ds = area)
    # we don't actually need 'tile' in mfhGrid; this is only for easier identification when debugging
    mmhGrid <- model.matrix(endemic, mfhGrid)
    # exclude intercept from endemic model matrix below, will be treated separately
    if (nbeta0 > 0) mmhGrid <- mmhGrid[,-1,drop=FALSE]

    # Extract endemic model components
    offsetGrid <- model.offset(mfhGrid)
    gridBlocks <- mfhGrid[["(BLOCK)"]]
    ds <- mfhGrid[["(ds)"]]


    ### Parse interaction functions

    if (hase) {

        ## Check interaction functions
        siaf <- do.call(".parseiaf", args = alist(siaf, "siaf", verbose=trace>0))
        constantsiaf <- attr(siaf, "constant")
        if (siaf$npars != nsiafpars) {
            stop("length of 'siafpars' (", nsiafpars,
                 ") does not match the 'siaf' specification (", siaf$npars, ")")
        }
        
        tiaf <- do.call(".parseiaf", args = alist(tiaf, "tiaf", verbose=trace>0))
        constanttiaf <- attr(tiaf, "constant")
        if (constanttiaf) gmax <- 1L
        if (tiaf$npars != ntiafpars) {
            stop("length of 'tiafpars' (", ntiafpars,
                 ") does not match the 'tiaf' specification (", tiaf$npars, ")")
        }

        ## Check control.siaf
        if (constantsiaf) control.siaf <- NULL else {
            stopifnot(is.null(control.siaf) || is.list(control.siaf))
        }

        ## Define function that integrates the two-dimensional 'siaf' function
        ## over the influence regions of the events
        if (!constantsiaf && !is.null(siaf$Fcircle) && !is.null(siaf$effRange))
        {
            ## pre-compute effective range of the 'siaf' (USED BY .siafInt)
            effRangeTypes <- rep_len(siaf$effRange(siafpars), nTypes)
        }
        .siafInt <- .siafIntFUN(siaf = siaf, noCircularIR = FALSE)
                                             # not certain beforehand
        .siafInt.args <- c(list(siafpars), control.siaf$F)

        ## Check gmax
        if (is.null(gmax)) {
            gmax <- max(tiaf$g(rep.int(0,nTypes), tiafpars, 1:nTypes))
            cat("assuming gmax =", gmax, "\n")
        } else if (!isScalar(gmax)) {
            stop("'gmax' must be scalar")
        }

    } else {
        if (!missing(siaf) && !is.null(siaf))
            warning("'siaf' can only be modelled in conjunction with an 'epidemic' process")
        if (!missing(tiaf) && !is.null(tiaf))
            warning("'tiaf' can only be modelled in conjunction with an 'epidemic' process")
        siaf <- tiaf <- NULL
        control.siaf <- NULL
    }


    ### print some information on the upcoming simulation

    txtPrehistory <- if (Nout == 0L) "no prehistory" else paste(Nout,
        ngettext(Nout, "event", "events"), "in the prehistory")
    cat("\nSimulating a", if (length(unpredMarks) > 2L) "marked",
        "spatio-temporal point pattern with",
        "\n\t-", nTypes, ngettext(nTypes, "event type", "event types"),
        "\n\t-", txtPrehistory)

    coefs <- c(
        if (nbeta0 > 1L) {
            setNames(beta0, paste0("h.type",typeNames))
        } else if (nbeta0 == 1L) setNames(beta0, "h.(Intercept)"),
        if (p > 0L) setNames(beta, paste("h",colnames(mmh),sep=".")),
        if (hase) setNames(gamma, paste("e",colnames(mme),sep=".")),
        if (hassiafpars) setNames(siafpars, paste("e.siaf",1:nsiafpars,sep=".")),
        if (hastiafpars) setNames(tiafpars, paste("e.tiaf",1:ntiafpars,sep="."))
    )
    cat("\n\t-", length(coefs), "coefficients:\n\n")
    print(coefs)



    ##########################################
    ### CIF of the temporal ground process ###
    ##########################################


    ### calculate integral of endemic component over W (= union of tiles)
    ### and over types for all time blocks in stgrid

    hIntWK <- if (hash) {
        dsexpeta <- local({
            eta <- drop(mmhGrid %*% beta)  # =0 if p = 0
            if (!is.null(offsetGrid)) eta <- offsetGrid + eta
            ds * exp(unname(eta))
        })
        fact <- if (nbeta0 > 1L) sum(exp(beta0)) else if (nbeta0 == 1L) nTypes*exp(unname(beta0)) else nTypes
        fact * c(tapply(dsexpeta, gridBlocks, sum))
    } else setNames(numeric(nrow(blockstarts)), blockstarts[,1]) # zeroes
    #<- is a named vector with names referencing BLOCK in stgrid


    ### helper function evaluating the epidemic terms of the ground intensity
    ### for a specific set of events (the lambdag function uses eTerms)
    
    eTermsCalc <- function (eventData, eventCoords)
    {
        # extract some marks from the eventData (USED INSIDE .siafInt() BELOW!)
        eventTypes <- as.integer(eventData$type)
        eps.s <- eventData$eps.s
        # distance to the border (required for siafInt below, and for epidataCS)
        bdist <- bdist(eventCoords, Wedges)
        # spatial influence regions of the events
        influenceRegion <- if (nrow(eventCoords) > 0L) .influenceRegions(
            events = SpatialPointsDataFrame(
                coords = eventCoords,
                data = data.frame(eps.s = eps.s, .bdist = bdist),
                match.ID = FALSE
            ),
            W = Wowin, npoly = nCircle2Poly, maxExtent = maxExtentOfW,
            clipper = "polyclip"
        ) else list()
        # epidemic terms
        if (!hase) {
            return(list(matrix(NA_real_, length(influenceRegion), 3L),
                        bdist, influenceRegion))
        }
        # epidemic model matrix (will be multiplied with gamma)
        mme <- buildmme(eventData)
        # integrate the two-dimensional 'siaf' function over the influence region
        siafInts <- if (length(influenceRegion) == 0L) numeric(0L) else {
            environment(.siafInt) <- environment()
            do.call(".siafInt", .siafInt.args)
        }
        # Matrix of terms in the epidemic component
        eTerms <- cbind(
            qSum = qSumTypes[eventTypes],
            expeta = epilinkinv(drop(mme %*% gamma)),
            siafInt = siafInts
        )
        # Return
        list(eTerms, bdist, influenceRegion)
    }


    ### function calculating the (upper bound) intensity of the ground process
    ### it relies on several objects for the epidemic component which are updated alongside simulation

    # t will be one of the break points in stgrid or an event time
    lambdagVec <- function (t, upper=FALSE)
    {
        ## endemic part
        hIntWKt <- hIntWK[[as.character(tBLOCK)]]

        ## epidemic part
        ejIntWt <- if (!hase || length(infectives) == 0L) numeric(0L) else {
            eTerms <- eTerms[infectives,,drop=FALSE]
            gTerm <- if (upper) {
                    rep.int(gmax, length(infectives))
                } else {
                    times <- eventMatrix[infectives,"time"]
                    types <- eventMatrix[infectives,"type"]
                    tiaf$g(t-times, tiafpars, types)
                }
            # ejIntWt only for infectives, others have 0
            setNames(apply(cbind(eTerms,gTerm), 1, prod), infectives)
        }

        c("0"=hIntWKt, ejIntWt)   # endemic component has index "0" !
    }


    ### helper function calculating the integral of lambdag from oldct to ct
    ### during simulation; it depends on the current values of the simulation

    add2Lambdag <- if (!hase || constanttiaf) {
            function () lambdagUpper * (ct-oldct)
        } else function () {
            # old endemic ground intensity * passed time
            hIntWKInt_oldct_ct <- lambdaghe[1L] * (ct-oldct)
            # integrated epidemic ground intensities of infectives (from oldct)
            ejIntWInt_oldct_ct <- if (length(infectives) == 0L) numeric(0L) else {
                eTermsProd <- apply(eTerms[infectives,,drop=FALSE], 1, prod)
                # integral of \id_{(0;eps.t]}(t-t_j) g(t-t_j \vert \kappa_j) from oldct to ct, for j in infectives
                # we can ignore the indicator because t-t_j is not >eps.t if t in [oldct;ct], because recoveries are change points
                times <- eventMatrix[infectives,"time"]
                types <- eventMatrix[infectives,"type"]
                gInt_0_ct    <- tiaf$G(ct   -times, tiafpars, types)
                gInt_0_oldct <- tiaf$G(oldct-times, tiafpars, types)
                gInt_oldct_ct <- gInt_0_ct - gInt_0_oldct
                eTermsProd * gInt_oldct_ct
            }
            sum(hIntWKInt_oldct_ct, ejIntWInt_oldct_ct)
        }



    ##################
    ### Simulation ###
    ##################

    ### Initialise values for simulation loop

    # all necessary components for an epidataCS object will be build along the simulation
    # let's start with the events of the prehistory
    tmp <- eTermsCalc(eventData, eventCoords)
    eTerms <- tmp[[1]]; rownames(eTerms) <- NULL
    bdists <- tmp[[2]]
    influenceRegions <- tmp[[3]]
    sources <- rep.int(list(integer(0L)), Nout)

    # Transform eventData into a matrix, which is faster with rbind
    # (factors will be recreated at the end of simulation)
    # simulated events will be subsequently appended to this matrix
    eventMatrix <- if (Nout == 0L) {
        matrix(numeric(0L), nrow=0L, ncol=ncol(eventData), dimnames=list(NULL, names(eventData)))
    } else {
        sapply(eventData, as.numeric, simplify = TRUE)    # prehistory
    }
    if (Nout == 1L) eventMatrix <- t(eventMatrix)
    # we will also know about the source of infection and corresponding BLOCK in stgrid
    navec <- rep.int(NA_real_, Nout)
    eventMatrix <- cbind(eventMatrix, source = navec, lambda.h = navec,
                         lambda.e = navec, Lambdag = navec, BLOCK = navec)

    # row indices of currently infective individuals
    infectives <- seq_len(Nout)

    # maximum total number of events (including prehistory)
    maxEvents <- Nout + nEvents

    # change points of lambdag
    stgridbreaks <- blockstarts[-1,2]
    Rtimes <- setNames(eventMatrix[,"time"]+eventMatrix[,"eps.t"],
                       infectives) # name indexes row of eventMatrix

    # index of next event (row in eventMatrix)
    j <- Nout + 1L

    # allocation of large objects for faster filling-in of new events
    allocated <- Nout
    ncolEventMatrix <- ncol(eventMatrix)
    newAllocation <- expression({
        eventMatrix <- rbind(eventMatrix, matrix(NA_real_, nrow = .allocate, ncol = ncolEventMatrix))
        eventCoords <- rbind(eventCoords, matrix(NA_real_, nrow = .allocate, ncol = 2L))
        eTerms <- rbind(eTerms, matrix(NA_real_, nrow = .allocate, ncol = 3L))
        bdists <- c(bdists, rep.int(NA_real_,.allocate))
        influenceRegions <- c(influenceRegions, vector(.allocate, mode="list"))
        sources <- c(sources, vector(.allocate, mode="list"))
        allocated <- allocated + .allocate
    })

    # current time point
    ct <- t0

    # current value of the cumulative intensity function of the ground process
    Lambdag <- 0

    # last point rejected?
    pointRejected <- FALSE

    # did we have numerical problems simulating from Exp(lambdagUpper) in the current loop?
    hadNumericalProblems0 <- FALSE

    # index of the current loop
    loopCounter <- 0L


    ### Let's Rock 'n' Roll

    if (trace > 0L) {
        cat("\nSimulation path (starting from t=", t0, "):\n---\n", sep="")
    } else {
        cat("\nSimulating (starting from t=", t0, ") ...\n", sep="")
    }

    while(j <= maxEvents && ct < T && (hash || length(infectives) > 0L))
    {
        loopCounter <- loopCounter + 1L
        if (trace > 0L && loopCounter %% trace == 0L) {
            cat(loopCounter, "@t =", ct, ":\t#simulated events =", j-1L-Nout,
            "\t#currently infective =", length(infectives),
            if (hase && !constanttiaf) paste("\tlast rejected?", pointRejected), "\n")
            flush.console()   # affects Windows only
        }

        # check if we need to allocate larger matrices
        if (j > allocated) {
            eval(newAllocation)
        }

        if (!pointRejected)  # what we have to do in the usual case
        {
            # we need the time block of stgrid corresponding to the new covariates,
            # i.e. search BLOCK such that t in [start; stop)
            tBLOCK <- blockstarts[findInterval(ct, blockstarts[,2]), 1]

            # Compute new infection intensity (upper bound)
            lambdaghe <- lambdagVec(ct, upper=TRUE)
            lambdagUpper <- sum(lambdaghe)

            # Determine time of next external change point
            changePoints <- c(nextblock = if (length(stgridbreaks) > 0L) stgridbreaks[1L],
                              Rtimes)
            nextChangePoint <- if (length(changePoints) > 0L) {
                   changePoints[which.min(changePoints)]    # don't use min() because need names
                } else Inf
        }
        pointRejected <- FALSE

        ## Simulate waiting time for the subsequent infection
        if (is.na(lambdagUpper)) {
            warning("simulation stopped due to undefined intensity")
            break
        }
        if (lambdagUpper < 0) {
            warning("simulation stopped due to negative overall intensity")
            break
        }
        Delta <- if (lambdagUpper == 0) Inf else tryCatch(
            rexp(1, rate = lambdagUpper),
            warning = function (w) { # rate was too small (for R >= 2.7.0,
                                     # rexp(1, Inf) returns 0 without warning)
                assign("hadNumericalProblems0", TRUE, inherits = TRUE)
                Inf
            })

        # Stop if lambdaStarMax too big meaning Delta == 0 (=> concurrent events)
        if (Delta == 0) {
            warning("simulation stopped due to infinite overall intensity")
            break
        }
        # Stop at all costs if end of simulation time [t0; T) has been reached
        if (isTRUE(min(ct+Delta, nextChangePoint) >= T)) {
            # ">=" because we don't want an event at "end"
            break
        }

        oldct <- ct
        if (ct + Delta > nextChangePoint) {
        ## Simulated time point is beyond the next time of intensity change (removal or endemic covariates)
            ct <- unname(nextChangePoint)
            # update cumulative intensity of the ground processes up to time ct,
            # i.e. add integral of lambdag from oldct to ct
            Lambdag <- Lambdag + add2Lambdag()
            # is this change point due to next time block in stgrid?
            if (names(nextChangePoint) == "nextblock") {
                stgridbreaks <- stgridbreaks[-1]
            } else { # i.e. change point due to recovery
                recoverer <- names(nextChangePoint)
                # update set of infectives
                infectives <- setdiff(infectives, recoverer)
                # remove recovery time from Rtimes
                .Rtimesidx <- match(recoverer, names(Rtimes))
                Rtimes <- Rtimes[-.Rtimesidx]
            }
        } else {
        ## Simulated time point lies within the thinning period
            ct <- ct + Delta
            # rejection sampling if non-constant temporal interaction kernel g
            if (hase && !constanttiaf) {
                # Calculate actual ground intensity for rejection probability at new ct
                lambdaghe <- lambdagVec(ct, upper=FALSE)
                lambdag <- sum(lambdaghe)
                # rejection sampling step
                if (lambdag/lambdagUpper < runif(1)) {
                    pointRejected <- TRUE
                    next
                }
            }
            # At this point, we have an actual event!

            # update cumulative intensity of the ground processes up to time ct,
            # i.e. add integral of lambdag from oldct to ct
            Lambdag <- Lambdag + add2Lambdag()
            # note that lambdaghe[1L] did not change by the above update in case of !constanttiaf,
            # which is expected by add2Lambdag (which requires the value of lambdag.h(oldct))

            # Where did the event come from: imported case or infection?
            .eventSource <- as.integer(sample(names(lambdaghe), 1L, prob=lambdaghe))

            # We now sample type and location
            if (.eventSource == 0L) {   # i.e. endemic source of infection
                .eventType <- sample(typeNames, 1L,
                                     prob=if (nbeta0 > 1L) exp(beta0))
                stgrididx <- which(gridBlocks == tBLOCK)
                .eventTile <- sample(stgrid$tile[stgrididx], 1L,
                                     prob=dsexpeta[stgrididx])  # this is a factor
                ## spsample doesn't guarantee that the sample will consist of
                ## exactly n points. if no point is sampled (very unlikely
                ## though), there would be an error
                ntries <- 1L
                .nsample <- 1L
                while(
                inherits(eventLocationSP <- try(
                    spsample(tiles[as.character(.eventTile),],
                             n=.nsample, type="random"),
                    silent = TRUE), "try-error")) {
                    .nsample <- 10L   # this also circumvents a bug in sp 1.0-0
                                      # (missing drop=FALSE in sample.Spatial())
                    if (ntries >= 1000) {
                        stop("'sp::spsample()' didn't succeed in sampling a ",
                             "point from tile \"", as.character(.eventTile), "\"")
                    }
                    ntries <- ntries + 1L
                }
                .eventLocation <- coordinates(eventLocationSP)[1L,,drop=FALSE]
            } else {    # i.e. source is one of the currently infective individuals
                sourceType <- eventMatrix[.eventSource,"type"]
                sourceCoords <- eventCoords[.eventSource,,drop=FALSE]
                sourceIR <- influenceRegions[[.eventSource]]
                sourceEpss <- eventMatrix[.eventSource,"eps.s"]
                .upperRange <- min(sourceEpss, maxExtentOfW)
                .eventType <- sample(typeNames[qmatrix[sourceType,]], 1L)
                .eventTypeCode <- match(.eventType, typeNames)
                eventLocationIR <- if (constantsiaf) {
                    as.matrix(coords.ppp(runifpoint(1L, win=sourceIR)))
                } else {
                    eventInsideIR <- FALSE
                    ntries <- 0L
                    while(!eventInsideIR) {
                        if (ntries >= 1000) {
                            stop("event location sampled by siaf$simulate() was",
                                 " rejected 1000 times (not in influence region)")
                        }
                        ntries <- ntries + 1L
                        eventLocationIR <- siaf$simulate(1L, siafpars,
                                                         .eventTypeCode,
                                                         .upperRange)
                        eventInsideIR <- inside.owin(eventLocationIR[,1],
                                                     eventLocationIR[,2],
                                                     sourceIR)
                    }
                    eventLocationIR
                }
                .eventLocation <- sourceCoords + eventLocationIR
                whichTile <- over(SpatialPoints(.eventLocation,
                                                proj4string=tiles@proj4string),
                                  tiles)
                if (is.na(whichTile)) {
                    warning("event generated at (",
                            paste(.eventLocation, collapse=","),
                            ") not in 'tiles'")
                    stop("'tiles' must cover all of 'W'")
                }
                .eventTile <- row.names(tiles)[whichTile]
                .eventTile <- factor(.eventTile, levels=tileLevels)
                if (is.na(.eventTile))
                    stop("tile \"", row.names(tiles)[whichTile],
                         "\" of simulated event is no level of stgrid$tile",
                         "\n-> verify row.names(tiles)")
            }
            .eventType <- factor(.eventType, levels=typeNames)

            # sample marks at this time and location
            .eventMarks <- rmarks(ct, .eventLocation)

            # gather event information
            .eventData <- data.frame(time=ct, tile=.eventTile, type=.eventType,
                .eventMarks, check.rows = FALSE, check.names = FALSE)

            # determine potential sources of infection (for epidataCS and lambda)
            .sources <- infectives[eventMatrix[infectives,"type"] %in% which(qmatrix[,.eventType])]
            if (length(.sources) > 0L) {
                .sdiffs <- .eventLocation[rep.int(1L,length(.sources)),,drop=FALSE] - eventCoords[.sources,,drop=FALSE]
                .sources <- .sources[sqrt(.rowSums(.sdiffs^2, length(.sources), 2L)) <= eventMatrix[.sources,"eps.s"]]
            }

            # calculate actual intensity at this time, location and type
            .mmhEvent <- buildmmh(.eventData)
            .etaEvent <- .mmhEvent %*% beta
            if (!is.null(.offsetEvent <- attr(.mmhEvent, "offset")))
                .etaEvent <- .etaEvent + .offsetEvent
            if (nbeta0 == 1L) {
                .etaEvent <- .etaEvent + beta0
            } else if (nbeta0 > 1L) {
                .etaEvent <- .etaEvent + beta0[.eventType]
            }
            .lambdah <- exp(.etaEvent)
            .lambdae <- if (hase && length(.sources) > 0L) {
                .sdiffs <- .eventLocation[rep.int(1L,length(.sources)),,drop=FALSE] - eventCoords[.sources,,drop=FALSE]
                .fSources <- siaf$f(.sdiffs, siafpars, eventMatrix[.sources,"type"])
                .gSources <- tiaf$g(ct - eventMatrix[.sources,"time"], tiafpars, eventMatrix[.sources,"type"])
                sum(eTerms[.sources,"expeta"] * .fSources * .gSources)
            } else 0

            # calculate terms of the epidemic component e_j(t,s) of the new infective
            tmp <- eTermsCalc(.eventData, .eventLocation)

            # Update objects
            eventMatrix[j,] <- c(ct, as.numeric(.eventTile), as.numeric(.eventType),
                                 sapply(.eventMarks, as.numeric), .eventSource,
                                 .lambdah, .lambdae, Lambdag, tBLOCK)
            eventCoords[j,] <- .eventLocation
            eTerms[j,] <- tmp[[1]]
            bdists[j] <- tmp[[2]]
            influenceRegions[[j]] <- tmp[[3]][[1]]
            sources[[j]] <- .sources

            # Update set of infectives and recovery times
            infectives <- c(infectives, j)
            Rtimes <- c(Rtimes, setNames(ct + .eventMarks[["eps.t"]], j))

            # Increment next event iterator
            j <- j + 1L
        }
    }
    if (trace > 0L) cat("---\n")

    
    ### update T if simulation ended preterm

    if (j > maxEvents || (!hash && length(infectives) == 0L)) {
        T <- ct
        # clip stgrid to effective time range of simulation
        stgrid <- subset(stgrid, start <= T)
        if (j > maxEvents) {
            cat("Maximum number of events (nEvents=", nEvents,
                ") reached @t = ", T, "\n", sep="")
        } else { # epidemic-only model
            cat("Simulation has ended preterm (no more infectives)",
                "@t =", T, "with", j-1L-Nout, "simulated events.\n")
        }
    } else { # ct >= T or ct+Delta >= T
        cat("Simulation has ended @t =", T, "with", j-1L-Nout,
            "simulated events.\n")
    }



    ##############
    ### Return ###
    ##############


    ### Throw warning in case of numerical difficulties

    if (hadNumericalProblems0) {
        warning("occasionally, the overall infection rate was numerically equal to 0")
    }


    ### throw an error if no events have been simulated
    ## because SpatialPoints[DataFrame]() does not allow the empty set, try:
    ## SpatialPoints(coords = matrix(numeric(0), 0, 2), bbox=bbox(W))
    
    if (j-1L == Nout) {
        stop("no events have been simulated")
    }

    
    ### transform eventMatrix back into a data.frame with original factor variables

    cat("\nPreparing simulated events for \"epidataCS\" ...\n")
    preEventData <- eventData

    # drop unused entries (due to large pre-allocation) from objects
    seqAlongEvents <- seq_len(j-1L)
    eventData <- as.data.frame(eventMatrix[seqAlongEvents,,drop=FALSE])

    # rebuild factor variables
    for (idx in which(sapply(preEventData, is.factor))) {
        origlevels <- levels(preEventData[[idx]])
        eventData[[idx]] <- factor(eventData[[idx]], levels=seq_along(origlevels), labels=origlevels)
    }

    # transform integer columns to integer
    eventData[c("source","BLOCK")] <- lapply(eventData[c("source","BLOCK")], as.integer)


    ### Append additional columns for an epidataCS object

    # add endemic covariates at events
    stgrididx <- apply(eventData[c("BLOCK","tile")], 1, function (x) {
        ret <- with(stgrid, which(BLOCK==as.integer(x[1L]) & tile==x[2L]))
        if (length(ret) == 0L) NA_integer_ else ret
        #<- events of the prehistory have missing BLOCKs, thus return NA
    })
    stgridIgnoreCols <- match(c("BLOCK",
                                setdiff(obligColsNames_stgrid, "start")),
                              names(stgrid))
    eventData <- cbind(eventData, stgrid[stgrididx, -stgridIgnoreCols, drop = FALSE])
    rownames(eventData) <- seqAlongEvents

    # add hidden columns
    eventData$.obsInfLength <- with(eventData, pmin(T-time, eps.t))
    eventData$.sources <- sources[seqAlongEvents]
    eventData$.bdist <- bdists[seqAlongEvents]
    eventData$.influenceRegion <- influenceRegions[seqAlongEvents]
    attr(eventData$.influenceRegion, "nCircle2Poly") <- nCircle2Poly
    attr(eventData$.influenceRegion, "clipper") <- "polyclip"


    ### Construct "epidataCS" object

    events <- SpatialPointsDataFrame(
        coords = eventCoords[seqAlongEvents,,drop=FALSE], data = eventData,
        proj4string = W@proj4string, match.ID = FALSE
        #, bbox = bbox(W))   # the bbox of SpatialPoints is defined as the actual
                             # bbox of the points and is also updated every time 
                             # when subsetting the SpatialPoints object
                             # -> useless to specify it as the bbox of W
    )

    if (.onlyEvents) {
        cat("Done.\n")
        attr(events, "timeRange") <- c(t0, T)
        attr(events, "runtime") <- proc.time()[[3]] - ptm
        return(events)
    }

    epi <- list(events=events, stgrid=stgrid, W=W, qmatrix=qmatrix)


    ### Return object of class "simEpidataCS"

    cat("Done.\n")
    # append configuration of the model
    epi$bbox <- bbox(W)
    epi$timeRange <- c(t0, T)
    epi$formula <- list(
                   endemic = if (typeSpecificEndemicIntercept) {
                       update(formula(endemic), ~ (1|type) + .)   # re-add to the formula
                   } else formula(endemic),
                   epidemic = formula(epidemic),
                   siaf = siaf, tiaf = tiaf
                   )
    if (epilink != "log") # set as attribute only if non-standard link function
        attr(epi$formula$epidemic, "link") <- epilink
    # coefficients as a numeric vector to be compatible with twinstim-methods
    epi$coefficients <- coefs  #list(beta0=beta0, beta=beta, gamma=gamma,
                               #     siafpars=siafpars, tiafpars=tiafpars)
    epi$npars <- c(nbeta0=nbeta0, p=p, q=q, nsiafpars=nsiafpars, ntiafpars=ntiafpars)
    epi$control.siaf <- control.siaf    # for R0.simEpidataCS
    epi$call <- cl
    epi$runtime <- proc.time()[[3]] - ptm
    class(epi) <- c("simEpidataCS", "epidataCS")
    return(epi)
}



#############################################################################
### much more efficient simulation for endemic-only models
### where intensities are piecewise constant and independent from the history
#############################################################################

## auxiliary function to calculate the endemic intensity by spatio-temporal cell
## from the model environment of a "twinstim" fit
.hGrid <- function (modelenv)
{
    .beta0 <- rep_len(if (modelenv$nbeta0==0L) 0 else modelenv$beta0,
                      modelenv$nTypes)
    hGrid <- sum(exp(.beta0)) * eval(modelenv$hGridExpr, envir = modelenv)
    blockstartstop <- modelenv$histIntervals[
        match(modelenv$gridBlocks, modelenv$histIntervals$BLOCK), ]
    data.frame(blockstartstop, tile = modelenv$gridTiles, hGrid = hGrid,
               hInt = hGrid * modelenv$ds * modelenv$dt,
               row.names = NULL, check.rows = FALSE, check.names = FALSE)
}

## simulate events from the endemic component of a "twinstim" fit
## this simulates pure (s,t,k) data with the only extra column being "tile"
simEndemicEvents <- function (object, tiles)
{
    ## check arguments
    if (is.null(modelenv <- environment(object)))
        stop("no model environment -- re-fit or update() with 'model=TRUE'")
    tileLevels <- levels(modelenv$gridTiles)
    tiles <- check_tiles(tiles, levels = tileLevels,
                         areas.stgrid = modelenv$ds[seq_along(tileLevels)],
                         keep.data = FALSE)
    
    ## calculate endemic intensity by spatio-temporal cell
    lambdaGrid <- .hGrid(modelenv)
    
    ## simulate number of events by cell
    nGrid <- rpois(n = nrow(lambdaGrid), lambda = lambdaGrid[["hInt"]])
    nTotal <- sum(nGrid)
    
    ## sample time points
    tps <- mapply(
        FUN = runif,
        n = nGrid, min = lambdaGrid[["start"]], max = lambdaGrid[["stop"]],
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
    
    ## sample types
    beta0 <- coeflist.default(coef(object), object$npars)[["nbeta0"]]
    nTypes <- nrow(object$qmatrix)
    types <- if (nTypes == 1L) {
        rep.int(1L, nTotal)
    } else {
        sample.int(n = nTypes, size = nTotal, replace = TRUE,
                   prob = if (length(beta0) > 1L) exp(beta0))
    }
    
    ## put event times, tiles, and types in a data frame
    events <- data.frame(
        ##lambdaGrid[rep.int(seq_len(nrow(lambdaGrid)), nGrid), c("tile", "BLOCK")],
        time = unlist(tps, recursive = FALSE, use.names = FALSE),
        tile = rep.int(lambdaGrid[["tile"]], nGrid),
        type = factor(types, levels = seq_len(nTypes), labels = rownames(object$qmatrix)),
        row.names = NULL, check.rows = FALSE, check.names = FALSE
    )
    
    ## sample coordinates from tiles
    nByTile <- tapply(X = nGrid, INDEX = lambdaGrid["tile"], FUN = sum)
    xyByTile <- sapply(
        X = names(nByTile),
        FUN = function (tile) {
            n <- nByTile[tile]
            if (n > 0L)
                coordinates(spsample(x = tiles[tile,], n = n, type = "random", iter = 10))
            ## else NULL
        },
        simplify = FALSE, USE.NAMES = TRUE
    )
    
    ## set coordinates of events
    events <- SpatialPointsDataFrame(
        coords = do.call("rbind", xyByTile),
        data = events[order(events$tile),],
        proj4string = tiles@proj4string,
        match.ID = FALSE)

    ## order by time
    events <- events[order(events$time),]
    row.names(events) <- seq_along(events)
    events
}


####################################################
### some twinstim-methods for "simEpidataCS" objects
####################################################

### wrapper for R0.twinstim

R0.simEpidataCS <- function (object, trimmed = TRUE, ...)
{
    R0.twinstim(object, newevents=object$events@data, trimmed = trimmed, ...)
}


### wrapper for intensityplot.twinstim

as.twinstim.simEpidataCS <- function (x)
{
    m <- do.call("twinstim", c(
        formula(x),
        list(data = quote(x), control.siaf = x$control.siaf,
             optim.args = list(par=coef(x), fixed=TRUE),
             model = TRUE, cumCIF = FALSE, verbose = FALSE)
        ))
    components2copy <- setdiff(names(m), names(x))
    for (comp in components2copy) x[[comp]] <- m[[comp]]
    environment(x) <- environment(m)
    class(x) <- c("simEpidataCS", "epidataCS", "twinstim")
    x
}

intensityplot.simEpidataCS <- function (x, ...)
{
    if (is.null(environment(x))) {
        objname <- deparse(substitute(x))
        message("Setting up the model environment ...")
        x <- as.twinstim.simEpidataCS(x)
        try({
            assign(objname, x, envir=parent.frame())
            message("Note: added model environment to '", objname,
                    "' for future use.")
        }, silent=TRUE)
    }
    intensityplot.twinstim(x, ...)
}


### the residual process Lambda_g(t) is stored with the simulated events

residuals.simEpidataCS <- function (object, ...)
{
    setNames(object$events$Lambdag,
             row.names(object$events))[!is.na(object$events$Lambdag)]
}



################################################################################
# A 'simulate' method for objects of class "twinstim".
################################################################################

### FIXME: actually stgrid's of simulations might have different time ranges
###        when nEvents is active -> atm, simplify ignores this

.rmarks <- function (data, t0, T)
{
    observedMarks <- subset(marks.epidataCS(data, coords = FALSE),
                            subset = time > t0 & time <= T)
    if (nrow(observedMarks) == 0L) {
        message("Note: 'data' does not contain any events during ('t0';'T'],\n",
                "      'rmarks' thus samples marks from all of 'data$events'")
        observedMarks <- marks.epidataCS(data, coords = FALSE)
    }
    observedMarks <- observedMarks[match("eps.t", names(observedMarks)):ncol(observedMarks)]
    rm(list = "data", inherits = FALSE)  # to save memory (environment is kept)
    function (t, s, n = 1L) {
        as.data.frame(lapply(observedMarks, function (x)
            sample(na.omit(x), size = n, replace = TRUE)),
                      optional = TRUE)
    }
}

simulate.twinstim <- function (object, nsim = 1, seed = NULL, data, tiles,
    newcoef = NULL, rmarks = NULL, t0 = NULL, T = NULL, nEvents = 1e5,
    control.siaf = object$control.siaf,
    W = data$W, trace = FALSE, nCircle2Poly = NULL, gmax = NULL,
    .allocate = 500, simplify = TRUE, ...)
{
    ptm <- proc.time()[[3]]
    cl <- match.call()


    ### Determine seed (this part is copied from stats:::simulate.lm with
    ### Copyright (C) 1995-2012 The R Core Team)

    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }


    ### Few checks

    stopifnot(inherits(object, "twinstim"), inherits(data, "epidataCS"))
    stopifnot(isScalar(nsim), nsim > 0)
    nsim <- as.integer(nsim)
    if (is.null(t0)) t0 <- object$timeRange[1]
    if (is.null(T))  T  <- object$timeRange[2]
    if (is.null(nCircle2Poly))
        nCircle2Poly <- attr(data$events$.influenceRegion, "nCircle2Poly")


    ### Retrieve arguments for simulation

    endemic  <- formula(object)$endemic
    epidemic <- formula(object)$epidemic
    # we don't need any reference to the original formula environment
    environment(endemic) <- environment(epidemic) <- .GlobalEnv
    if (is.null(rmarks))
        rmarks <- .rmarks(data, t0 = t0, T = T)
    theta <- coef(object)
    if (!is.null(newcoef)) {
        newcoef <- check_twinstim_start(newcoef)
        newcoef <- newcoef[names(newcoef) %in% names(theta)]
        theta[names(newcoef)] <- newcoef
    }
    thetalist <- coeflist.default(theta, object$npars)


    ### Run the simulation(s)

    # establish call
    simcall <- call("simEpidataCS", endemic=endemic, epidemic=epidemic,
                    siaf=quote(formula(object)$siaf),
                    tiaf=quote(formula(object)$tiaf),
                    qmatrix=quote(object$qmatrix),
                    rmarks=quote(rmarks), events=quote(data$events),
                    stgrid=quote(data$stgrid), tiles=quote(tiles),
                    beta0=thetalist[[1L]], beta=thetalist[[2L]],
                    gamma=thetalist[[3L]], siafpars=thetalist[[4L]],
                    tiafpars=thetalist[[5L]], epilink = .epilink(object),
                    t0=t0, T=T, nEvents=nEvents,
                    control.siaf=control.siaf,
                    W=quote(W), trace=trace, nCircle2Poly=nCircle2Poly,
                    gmax=gmax, .allocate=.allocate,
                    .skipChecks=TRUE, .onlyEvents=FALSE)
    
    # First simulation
    if (nsim > 1L) {
        cat("\nTime at beginning of simulation:", as.character(Sys.time()), "\n")
        cat("Simulation 1 /", nsim, "...\n")
        cat("-------------------------------------------------------------------------------\n")
    }
    res <- eval(simcall)
    if (nsim > 1L) {
        cat("\n-------------------------------------------------------------------------------\n")
        cat("Runtime of first simulation:", res$runtime, "seconds\n")
        cat("Estimated finishing time:", as.character(Sys.time() + (nsim-1) * res$runtime), "\n\n")
        # set up list of simulations
        res <- if (simplify) {
                with(res, list(
                    eventsList=c(structure(events, timeRange = timeRange, runtime = runtime),
                                 vector(nsim-1L, mode="list")),
                    stgrid=stgrid, W=W, qmatrix=qmatrix, formula=formula,
                    coefficients=coefficients, npars=npars, call=call
                ))
            } else {
                c(list(res), vector(nsim-1L, mode="list"))
            }
        # force garbage collection
        gc()
        # run the remaining simulations
        simcall$.onlyEvents <- simplify
        for (i in 2:nsim) {
            cat("Simulation", sprintf(paste0("%",nchar(nsim),"i"), i), "/", nsim, "...")
            capture.output(resi <- eval(simcall))
            .nEvents <- if (simplify) sum(!is.na(resi$source)) else {
                sum(!is.na(resi$events$source))
            }
            .T <- if (simplify) attr(resi,"timeRange")[2] else resi$timeRange[2]
            cat("\tsimulated", .nEvents, "events", if (nEvents == .nEvents)
                "(reached maximum)", "up to time", .T, "\n")
            if (simplify) res$eventsList[[i]] <- resi else res[[i]] <- resi
        }
        cat("\nDone (", as.character(Sys.time()), ").\n", sep="")
    }
    attr(res, "call") <- cl
    attr(res, "seed") <- RNGstate
    attr(res, "runtime") <- proc.time()[[3]] - ptm
    class(res) <- if (nsim == 1L) {
            c("simEpidataCS", "epidataCS")
        } else {
            attr(res, "simplified") <- simplify
            c("simEpidataCSlist")
        }
    res
}



### print method for lists of simulated epidemics

print.simEpidataCSlist <- function (x, ...)
{
    cat("\nCall:\n")
    print.default(attr(x, "call"))
    simplified <- attr(x, "simplified")
    nsim <- if (simplified) length(x$eventsList) else length(x)
    cat("\n")
    cat(if (simplified) "Simplified list" else "List", "of", nsim,
        "simulated epidemics of class \"simEpidataCS\" (not printed)\n\n")
    invisible(x)
}

"[[.simEpidataCSlist" <- function (x, i) {
    simplified <- attr(x, "simplified")
    if (simplified) {
        x <- unclass(x)
        x$eventsList <- x$eventsList[[i]]
        names(x)[names(x) == "eventsList"] <- "events"
        x <- append(x, list(timeRange = attr(x$events, "timeRange")), after=4L)
        x$runtime <- attr(x$events, "runtime")
        attr(x$events, "timeRange") <- attr(x$events, "runtime") <- NULL
        class(x) <- c("simEpidataCS", "epidataCS")
        x
    } else NextMethod("[[")
}

plot.simEpidataCSlist <- function (x,
    which = NULL, mfrow = n2mfrow(length(which)),
    main = paste("Simulated epidemic", which),
    aggregate = c("time", "space"), subset, ...)
{
    simplified <- attr(x, "simplified")
    nsim <- if (simplified) length(x$eventsList) else length(x)
    if (is.null(which)) {
        which <- seq_len(nsim)
        if (nsim > 4) which <- sample(which, 4L)
    }
    opar <- par(mfrow = mfrow); on.exit(par(opar))
    main <- rep_len(main, length(which))
    for (i in seq_along(which)) {
        do.call("plot", args=list(x=quote(x[[which[i]]]), aggregate=aggregate,
                        subset=substitute(subset), main = main[i], ...))
    }
}
