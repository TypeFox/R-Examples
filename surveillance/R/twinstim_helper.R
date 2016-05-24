################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Some internal helper functions for "twinstim".
###
### Copyright (C) 2009-2015 Sebastian Meyer
### $Revision: 1285 $
### $Date: 2015-03-24 15:26:51 +0100 (Die, 24. Mär 2015) $
################################################################################


### Determines indexes of potential sources of infection

## determine potential sources of the i'th event
## all arguments but i and qmatrix are nEvents-vectors
## -> determine potential sources for eventTimes[i], eventsTypes[i] with
## distances distvec_j = ||s_i - s_j||
determineSources <- function (i, eventTimes, removalTimes, distvec, eps.s,
    eventTypes = NULL, qmatrix)
{
    tp <- eventTimes[i]
    infectivity <- (eventTimes < tp) & (removalTimes >= tp)
    #<- eventTimes<tp, not "=" because CIF is left-continuous.
    #Also guarantees no self-infection
    proximity <- distvec <= eps.s
    sources <- if (is.null(eventTypes)) {
        which(infectivity & proximity)
    } else {
        type <- eventTypes[i]
        typeInfective <- qmatrix[,type] # indexed by integer code of factor
        #<- logical vector indicating for each type if it could infect type of i
        matchType <- typeInfective[eventTypes]
        which(infectivity & proximity & matchType)
    }
    unname(sources)
}

## determine the .sources for an epidataCS object, i.e.
## lapply the previous function to all of object$events
determineSources.epidataCS <- function (object)
{
    eventTimes <- object$events@data$time
    removalTimes <- eventTimes + object$events@data$eps.t
    eventDists <- as.matrix(dist(object$events@coords, method = "euclidean"))
    lapply(seq_along(eventTimes), function (i) {
        determineSources(i, eventTimes, removalTimes, eventDists[i,],
                         object$events@data$eps.s, object$events@data$type,
                         object$qmatrix) 
    })
}



### Check matrix Q

checkQ <- function (qmatrix, typeNames)
{
    nTypes <- length(typeNames)
    qmatrix <- as.matrix(qmatrix)
    stopifnot(nrow(qmatrix) == ncol(qmatrix))
    if (is.null(dimnames(qmatrix))) {
        if (nrow(qmatrix) != nTypes) {
            stop("'qmatrix' must be a ", nTypes, "x", nTypes, " matrix")
        }
        dimnames(qmatrix) <- list(typeNames, typeNames)
    } else {
        stopifnot(rownames(qmatrix) == colnames(qmatrix))
        typesIdx <- match(typeNames, rownames(qmatrix), nomatch=NA_integer_)
        if (idx <- match(TRUE, is.na(typesIdx), nomatch=0L)) {
            stop("missing entries for type '", typeNames[idx], "' in 'qmatrix'")
        }
        qmatrix <- qmatrix[typesIdx,typesIdx,drop=FALSE]
    }
    storage.mode(qmatrix) <- "logical"   # convert entries to TRUE/FALSE by as.logical
    qmatrix
}


### Get row index of 'stgrid' where an event is located (spatio-temporally)
### Here, search BLOCK such that t in (start;stop], i.e. an event at 'stop' is
### still attributed to the previous interval

gridcellOfEvent <- function (t, tilename, stgrid)
{
    ## idx <- with(stgrid, which(tile == tilename & start < t & stop >= t))
    
    ## ~5x faster alternative assuming a full BLOCK x tile grid, which is
    ## sorted by BLOCK and tile (tile varying first), specifically there must be
    ## all levels(stgrid$tile) in every BLOCK in that order;
    ## this structure is guaranteed by check_stgrid()
    blockstart <- match(TRUE, stgrid$stop >= t)
    idx <- blockstart + match(tilename, levels(stgrid$tile)) - 1L
    
    lidx <- length(idx)
    if (lidx == 0L) NA_integer_ else if (lidx == 1L) idx else {
        stop("'stgrid' has overlapping spatio-temporal grid cells")
    }
}


## Crude estimate for a start value of the endemic intercept
## assuming the model only had a single-cell endemic component
## (rate of homogeneous Poisson process scaled for the offset)
crudebeta0 <- function (nEvents, offset.mean, W.area, period, nTypes)
{
    ## nEvents = exp(offset + beta0) * W.area * period * nTypes
    log(nEvents/W.area/period/nTypes) - offset.mean
}


### Really internal helper function, which constructs the function that
### integrates the two-dimensional 'siaf' function over the influence regions of
### the events. The only argument of the returned function is 'siafpars'.
### The returned function is defined in the callers environment, where the
### variables used in the function are available (inside twinstim() or
### simEpidataCS()).

.siafIntFUN <- function (siaf,
    noCircularIR, #= all(eps.s>bdist) = all(sapply(influenceRegion, function(x)
                  #                           is.null(attr(x,"radius"))))
    parallel = FALSE
){
    ## the following variables are unused here, because the environment of
    ## FUN will be set to the parent.frame(), where the variables exist
    ## they are only included to avoid the notes in R CMD check 
    iRareas <- influenceRegion <- eventTypes <- eps.s <- bdist <- effRanges <- NULL
    
    ## define the siaf integration function depending on the siaf specification 
    FUN <- if (attr(siaf, "constant"))
    {
        if (exists("iRareas", where=parent.frame(), mode="numeric")) {
            ## in twinstim(), 'iRareas' are pre-defined to save
            ## computation time (data are fixed during fitting)
            function (siafpars) iRareas
        } else {
            function (siafpars)
                vapply(X = influenceRegion, FUN = attr, which = "area",
                       FUN.VALUE = 0, USE.NAMES = FALSE)
        }
    } else if (is.null(siaf$Fcircle) || # if siaf$Fcircle not available
               (is.null(siaf$effRange) && noCircularIR))
    {
        ## Numerically integrate 'siaf' over each influence region
        mapplyFUN(
            c(alist(siaf$F, influenceRegion, type=eventTypes),
              list(MoreArgs=quote(list(siaf$f, siafpars, ...)),
                   SIMPLIFY=TRUE, USE.NAMES=FALSE)),
            ##<- we explicitly quote() the ...-part instead of simply including
            ##   it in the above alist() - only to make checkUsage() happy
            parallel = parallel)
    } else if (is.null(siaf$effRange)) # use Fcircle but only delta-trick
    {
        mapplyFUN(
            c(alist(function (iR, type, eps, bdisti, siafpars, ...)
                    if (eps <= bdisti) # influence region completely inside W
                        siaf$Fcircle(eps, siafpars, type)
                    else # numerically integrate over influence region
                        siaf$F(iR, siaf$f, siafpars, type, ...)
                    ,
                    influenceRegion, eventTypes, eps.s, bdist),
              list(MoreArgs=quote(list(siafpars, ...)),
                   SIMPLIFY=TRUE, USE.NAMES=FALSE)),
            parallel = parallel)
    } else { # fast Fcircle integration considering the delta-trick AND effRange
        .ret <- mapplyFUN(
            c(alist(function (iR, type, eps, bdisti, effRange, siafpars, ...)
                    if (eps <= bdisti) # influence region completely inside W
                        siaf$Fcircle(eps, siafpars, type)
                    else if (effRange <= bdisti) # effective region inside W
                        siaf$Fcircle(bdisti, siafpars, type)
                    else # numerically integrate over influence region
                        siaf$F(iR, siaf$f, siafpars, type, ...)
                    ,
                    influenceRegion, eventTypes, eps.s, bdist, effRanges),
              list(MoreArgs=quote(list(siafpars, ...)),
                   SIMPLIFY=TRUE, USE.NAMES=FALSE)),
            ## before: compute computationally effective range of the 'siaf'
            ## for the current 'siafpars' for each event (type):
            before = expression(
            effRangeTypes <- rep_len(siaf$effRange(siafpars), nTypes),
            effRanges <- effRangeTypes[eventTypes]   # N-vector
            ),
            parallel = parallel)        
        if (exists("effRangeTypes", where=parent.frame(), mode="numeric")) {
            ## in simEpidataCS effRangeTypes is pre-calculated outside siafInt to
            ## save computation time ('siafpars' is constant during simulation)
            body(.ret)[[grep("^effRangeTypes <-", body(.ret))]] <- NULL
        }
        .ret
    }
    
    ## set the environment of the siafInt function to the callers environment
    ## (i.e. inside twinstim() or simEpidataCS())
    ## where the variables used in the function are defined
    environment(FUN) <- parent.frame()
    FUN
}


### Helper function, which constructs the function that integrates the 'tiaf'.
### The returned function is defined in the callers environment, where the
### variables used in the function are available (inside twinstim() or
### simEpidataCS()).

.tiafIntFUN <- function ()
{
    ## the following variables are unused here, because the environment of
    ## FUN will be set to the parent.frame(), where the variables exist
    ## they are only included to avoid the notes in R CMD check 
    gIntLower <- gIntUpper <- eventTypes <- tiaf <- NULL
    
    ## from, to and type may be vectors of compatible lengths
    FUN <- function(tiafpars, from = gIntLower, to = gIntUpper,
                    type = eventTypes, G = tiaf$G)
    {
        tiafIntUpper <- G(to, tiafpars, type)
        tiafIntLower <- G(from, tiafpars, type)
        tiafIntUpper - tiafIntLower
    }
    
    ## set the environment of the tiafInt function to the callers environment
    ## (i.e. inside twinstim() or simEpidataCS())
    ## where the default argument values are defined
    environment(FUN) <- parent.frame()
    FUN
}


### rename control arguments with optim names to have names compatible with nlminb

control2nlminb <- function (control, defaults)
{
    renamelist <- cbind(optim  = c("maxit", "REPORT", "abstol", "reltol"),
                        nlminb = c("iter.max", "trace", "abs.tol", "rel.tol"))
    for (i in which(renamelist[,"optim"] %in% names(control))) {
        fromname <- renamelist[i, "optim"]
        toname <- renamelist[i, "nlminb"]
        if (is.null(control[[toname]])) {
            control[[toname]] <- control[[fromname]]
        }
        control[[fromname]] <- NULL
    }
    defaults[names(control)] <- control
    defaults
}


### Helper for iaf-checks:
### Checks if FUN has three arguments (s/t, pars, type) and
### eventually adds the last two

.checknargs3 <- function (FUN, name)
{
    FUN <- match.fun(FUN)
    NARGS <- length(formals(FUN))
    if (NARGS == 0L) {
        stop("the function '", name, "' must accept at least one argument")
    } else if (NARGS == 1L) {
        formals(FUN) <- c(formals(FUN), alist(pars=, types=))
    } else if (NARGS == 2L) {
        formals(FUN) <- c(formals(FUN), alist(types=))
    }
    FUN
}


### Internal wrapper used in twinstim() and simEpidataCS() to evaluate the siaf
### and tiaf arguments. If succesful, returns checked interaction function.

.parseiaf <- function (iaf, type, eps = NULL, verbose = TRUE)
{
    type <- match.arg(type, choices=c("siaf", "tiaf"), several.ok=FALSE)
    res <- if (missing(iaf) || is.null(iaf)) {
        if (verbose) {
            message("assuming constant ",
                    switch(type, siaf="spatial", tiaf="temporal"),
                    " interaction '", type, ".constant()'")
        }
        do.call(paste(type, "constant", sep="."), args=alist())
    } else if (is.list(iaf)) {
        ret <- do.call(type, args = iaf)
        ## keep special attributes
        attr(ret, "knots") <- attr(iaf, "knots")
        attr(ret, "maxRange") <- attr(iaf, "maxRange")
        attr(ret, "Boundary.knots") <- attr(iaf, "Boundary.knots")
        attr(ret, "constant") <- attr(iaf, "constant")
        ret
    } else if (is.vector(iaf, mode = "numeric")) {
        do.call(paste(type,"step",sep="."), args = list(knots = iaf))
    } else {
        stop("'", as.character(substitute(iaf)),
             "' must be NULL (or missing), a list (-> continuous ",
             "function), or numeric (-> knots of step function)")
    }
    ## indicate if this is a constant iaf
    attr(res, "constant") <- isTRUE(attr(res, "constant"))
    ## attach unique interaction ranges
    if (!is.null(eps)) {         # in simEpidataCS() eps is not known beforehand
        attr(res, "eps") <- sort(unique(eps))
    }
    return(res)
}


### Construct a call/function for mapply or parallel::mcmapply, respectively
## args: alist() of arguments for mapply()
## before,after: expressions to be prepended/appended to the function body,
##               where "res" will be the result of mapply()

mapplyCall <- function (args, cores = 1L)
{
    parallel <- is.name(cores) || cores > 1L
    mapplyFUN <- if (parallel) quote(parallel::mcmapply) else quote(mapply)
    parallelArgs <- list(mc.preschedule=TRUE, mc.cores=cores)
    as.call(c(mapplyFUN, args, if (parallel) parallelArgs))
}

mapplyFUN <- function (args, before = list(), after = list(), parallel = TRUE)
{
    FUN <- as.function(alist(siafpars=, ...=, NULL),
                       envir=parent.frame())
    body(FUN) <- mapplyCall(args, if (parallel) quote(cores) else 1L)
    if (length(after) + length(before) > 0) {
        body(FUN) <- as.call(c(
            list(as.name("{")),
            before,
            if (length(after)) call("<-", as.name("res"), body(FUN)) else body(FUN),
            after))
    }
    FUN
}


### parse the list or vector of start values

check_twinstim_start <- function (start)
{
    if (is.null(start)) {
        return(start)
    } else if (is.list(start)) { # convert allowed list specification to vector
        stopifnot(names(start) %in% c("endemic", "epidemic", "h", "e",
                                      "siaf", "tiaf", "e.siaf", "e.tiaf"))
        names(start)[names(start) == "endemic"] <- "h"
        names(start)[names(start) == "epidemic"] <- "e"
        names(start)[names(start) == "siaf"] <- "e.siaf"
        names(start)[names(start) == "tiaf"] <- "e.tiaf"
        start <- unlist(start, recursive=FALSE, use.names=TRUE)
    }
    if (!is.vector(start, mode="numeric") || is.null(names(start)))
        stop("parameter values must be named and numeric")
    return(start)
}
