################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Standard S3-methods for "epidataCS" objects, which represent
### CONTINUOUS SPATIO-temporal infectious disease case data
###
### Copyright (C) 2009-2015 Sebastian Meyer
### $Revision: 1285 $
### $Date: 2015-03-24 15:26:51 +0100 (Die, 24. Mär 2015) $
################################################################################


### Number of events

nobs.epidataCS <- function (object, ...) length(object$events)


### UPDATE eps.s, eps.t, qmatrix OR nCircle2Poly IN AN EXISTING epidataCS OBJECT

# all arguments but 'object' are optional, the ... argument is unused
update.epidataCS <- function (object, eps.t, eps.s, qmatrix, nCircle2Poly, ...)
{
    nEvents <- nobs(object)

    # Check and update eps.t
    if (!missing(eps.t)) {
        stopifnot(is.numeric(eps.t), eps.t > 0)
        object$events$eps.t <- eps.t
    }

    # Initialise indicator of which influenceRegions to update
    ir2update <- logical(nEvents)   # all FALSE

    # Check and update eps.s
    if (!missing(eps.s)) {
        stopifnot(is.numeric(eps.s), eps.s > 0)
        oldeps.s <- object$events$eps.s
        object$events$eps.s <- eps.s
        ir2update <- oldeps.s != object$events$eps.s
    }

    # Check nCircle2Poly
    nCircle2Poly <- if (missing(nCircle2Poly)) {
        attr(object$events$.influenceRegion, "nCircle2Poly")
    } else {
        stopifnot(isScalar(nCircle2Poly))
        ir2update <- rep.int(TRUE, nEvents)
        as.integer(nCircle2Poly)
    }

    # Update influenceRegions of events
    if (any(ir2update)) {
        clipper <- attr(object$events$.influenceRegion, "clipper")
        if (is.null(clipper))  # epidataCS < 1.8-1
            clipper <- "polyclip"
        object$events$.influenceRegion[ir2update] <-
            .influenceRegions(object$events[ir2update,], object$W, nCircle2Poly,
                              clipper = clipper)
        attr(object$events$.influenceRegion, "nCircle2Poly") <- nCircle2Poly
    }

    # Check qmatrix
    if (!missing(qmatrix)) object$qmatrix <- checkQ(qmatrix, levels(object$events$type))

    #hoehle @ 16 Apr 2011 - bug fix. .obsInfLength was not handled
    # Update length of infection time, i.e. length = min(T-time, eps.t)
    if (!missing(eps.t)) {
      timeRange <- with(object$stgrid, c(start[1], stop[length(stop)]))
      object$events$.obsInfLength <- with(object$events@data, pmin(timeRange[2]-time, eps.t))
    }

    # Update .sources
    if (!missing(eps.t) || !missing(eps.s) || !missing(qmatrix)) {
        object$events$.sources <- determineSources.epidataCS(object)
    }

    # Done update.
    return(object)
}



### subsetting epidataCS, i.e. select only part of the events,
### but retain stgrid and W. If any event types disappear due to subsetting,
### these types will be dropped from the factor levels and from qmatrix

"[.epidataCS" <- function (x, i, j, ..., drop = TRUE)
{
    ## rescue attributes of .influenceRegion (dropped when indexing)
    iRattr <- attributes(x$events$.influenceRegion)
    
    ## apply [,SpatialPointsDataFrame-method (where "drop" is ignored)
    cl <- sys.call()
    cl[[1]] <- as.name("[")
    cl[[2]] <- substitute(x$events)
    x$events <- eval(cl, envir=parent.frame())
    
    ## assure valid epidataCS after subsetting
    if (!missing(j)) {                # only epidemic covariates may be selected
        endemicVars <- setdiff(names(x$stgrid), c(
            reservedColsNames_stgrid, obligColsNames_stgrid))
        if (!all(c(reservedColsNames_events, obligColsNames_events, 
                   endemicVars) %in% names(x$events))) {
            stop("only epidemic covariates may be removed from 'events'")
        }
    }
    if (!missing(i)) {
        ## update .sources
        x$events$.sources <- determineSources.epidataCS(x)
        if (drop) {
            ## update type levels and qmatrix (a type could have disappeared)
            x$events$type <- x$events$type[drop=TRUE]
            typeNames <- levels(x$events$type)
            if (!identical(rownames(x$qmatrix), typeNames)) {
                message("Note: dropped type(s) ",
                        paste0("\"", setdiff(rownames(x$qmatrix), typeNames), "\"",
                               collapse = ", "))
                x$qmatrix <- checkQ(x$qmatrix, typeNames)
            }
        }
    }
    
    ## restore attributes of .influenceRegion
    attributes(x$events$.influenceRegion) <- iRattr
    
    ## done
    return(x)
}


## The subset method for epidataCS-objects is adapted from
## base::subset.data.frame (authored by Peter 
## Dalgaard and Brian Ripley, Copyright (C) 1995-2012
## The R Core Team) with slight modifications only
## (we just replace 'x' by 'x$events@data' for evaluation of subset and select)

subset.epidataCS <- function (x, subset, select, drop = TRUE, ...)
{
    if (missing(subset)) 
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, x$events@data, parent.frame()) # HERE IS A MOD
        if (!is.logical(r)) stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    if (missing(select)) 
        vars <- TRUE
    else {
        nl <- as.list(seq_along(x$events@data)) # HERE IS A MOD
        names(nl) <- names(x$events@data)       # HERE IS A MOD
        vars <- eval(substitute(select), nl, parent.frame())
    }
    x[r, vars, drop = drop]       # this calls the [.epidataCS-method from above
}


## Subset epidataCS object using head and tail methods (which use [.epidataCS)

head.epidataCS <- function (x, n = 6L, ...)
    head.matrix(x, n = n, ...)

tail.epidataCS <- function (x, n = 6L, ...)
{
    # ugly hack for tail.matrix because I don't want to register a
    # dim-method for class "epidataCS"
    nrow <- function (x) base::nrow(x$events)
    my.tail.matrix <- tail.matrix
    environment(my.tail.matrix) <- environment()
    ##<- such that the function uses my local nrow definition
    my.tail.matrix(x, n = n, addrownums=FALSE, ...)
}



### extract marks of the events (actually also including time and tile)

idxNonMarks <- function (x)
{
    endemicCovars <- setdiff(names(x$stgrid), c(
        reservedColsNames_stgrid, obligColsNames_stgrid))
    match(c(reservedColsNames_events, endemicCovars), names(x$events@data))
}

marks.epidataCS <- function (x, coords = TRUE, ...)
{
    if (coords) { # append coords (cp. as.data.frame.SpatialPointsDataFrame)
        data.frame(x$events@data[-idxNonMarks(x)], x$events@coords)
    } else { # return marks without coordinates
        x$events@data[-idxNonMarks(x)]
    }
}



### permute event times and/or locations holding remaining columns fixed

permute.epidataCS <- function (x, what = c("time", "space"), keep)
{
    stopifnot(inherits(x, "epidataCS"))
    what <- match.arg(what)
    
    ## permutation index
    perm <- if (missing(keep)) {
        sample.int(nobs.epidataCS(x))
    } else { # some events should not be relabeled
        keep <- eval(substitute(keep), envir = x$events@data,
                     enclos = parent.frame())
        stopifnot(is.logical(keep), !is.na(keep))
        which2permute <- which(!keep)
        howmany2permute <- length(which2permute)
        if (howmany2permute < 2L) {
            message("Note: data unchanged ('keep' all)")
            return(x)
        }
        perm <- seq_len(nobs.epidataCS(x))
        perm[which2permute] <- which2permute[sample.int(howmany2permute)]
        perm
    }
    
    ## rescue attributes of .influenceRegion (dropped when indexing)
    iRattr <- attributes(x$events@data$.influenceRegion)
    
    ## permute time points and/or locations
    PERMVARS <- if (what == "time") {
        c("time", "BLOCK", "start", ".obsInfLength")
    } else {
        x$events@coords <- x$events@coords[perm,,drop=FALSE]
        c("tile", ".bdist", ".influenceRegion")
    }
    x$events@data[PERMVARS] <- x$events@data[perm, PERMVARS]
    
    ## re-sort on time if necessary
    if (what == "time") {
        x$events <- x$events[order(x$events@data$time), ]
    }
    
    ## .sources and endemic variables need an update
    x$events@data$.sources <- determineSources.epidataCS(x)
    ENDVARS <- setdiff(names(x$stgrid),
                       c(reservedColsNames_stgrid, obligColsNames_stgrid))
    gridcellsOfEvents <- match(
        do.call("paste", c(x$events@data[c("BLOCK", "tile")], sep = "\r")),
        do.call("paste", c(x$stgrid[c("BLOCK", "tile")], sep = "\r"))
    )
    x$events@data[ENDVARS] <- x$stgrid[gridcellsOfEvents, ENDVARS]
    
    ## restore attributes of .influenceRegion
    attributes(x$events@data$.influenceRegion) <- iRattr
    
    ## done
    x
}



### printing methods

print.epidataCS <- function (x, n = 6L, digits = getOption("digits"), ...)
{
    print.epidataCS_header(
        timeRange = c(x$stgrid$start[1L], x$stgrid$stop[nrow(x$stgrid)]),
        bbox = bbox(x$W),
        nBlocks = length(unique(x$stgrid$BLOCK)),
        nTiles = nlevels(x$stgrid$tile),
        digits = digits
        )
    cat("Types of events: ")
    str(levels(x$events$type), give.attr = FALSE, give.head = FALSE,
        width = getOption("width") - 17L)
    cat("Overall number of events:", nEvents <- nobs(x), "\n\n")
    
    visibleCols <- grep("^\\..+", names(x$events@data), invert = TRUE)
    if (nEvents == 0L) { # not handled by [,SpatialPointsDataFrame-method
                         # and thus actually not supported by "epidataCS"
        ## display header only
        print(data.frame(coordinates = character(0L), x$events@data[visibleCols]))
    } else {
        ## 2014-03-24: since sp 1.0-15, print.SpatialPointsDataFrame()
        ## appropriately passes its "digits" argument to print.data.frame()
        print(head.matrix(x$events[visibleCols], n = n), digits = digits, ...)
        if (n < nEvents) cat("[....]\n")
    }
    invisible(x)
}

print.epidataCS_header <- function (timeRange, bbox, nBlocks, nTiles,
                                    digits = getOption("digits"))
{
    bboxtxt <- paste(
        apply(bbox, 1, function (int) paste0(
            "[", paste(format(int, trim=TRUE, digits=digits), collapse=", "), "]"
            )),
        collapse = " x ")
    cat("Observation period:",
        paste(format(timeRange, trim=TRUE, digits=digits),
              collapse = " - "), "\n")
    cat("Observation window (bounding box):", bboxtxt, "\n")
    cat("Spatio-temporal grid (not shown):",
        nBlocks, ngettext(nBlocks, "time block,", "time blocks"),
        "x",
        nTiles, ngettext(nTiles, "tile", "tiles"),
        "\n")
}


### SUMMARY
# the epidemic is summarized by the following returned components:
# timeRange, nEvents, eventTimes, eventCoords, nSources, as well as
# - tile/typetable: number of events per tile/type
# - counter: number of infective individuals as stepfun

summary.epidataCS <- function (object, ...)
{
    res <- list(
        timeRange = with(object$stgrid, c(start[1], stop[length(stop)])),
        bbox = bbox(object$W),
        nBlocks = length(unique(object$stgrid$BLOCK)),
        nEvents = nobs(object),
        nTypes = nlevels(object$events$type),
        eventTimes = object$events$time,
        eventCoords = coordinates(object$events),
        eventTypes = object$events$type,
        eventRanges = object$events@data[c("eps.t", "eps.s")],
        eventMarks = marks.epidataCS(object), 
        tileTable = c(table(object$events$tile)),
        typeTable = c(table(object$events$type)),
        counter = as.stepfun.epidataCS(object),
        nSources = sapply(object$events$.sources, length)
        )
    class(res) <- "summary.epidataCS"
    res
}

print.summary.epidataCS <- function (x, ...)
{
    print.epidataCS_header(timeRange = x$timeRange, bbox = x$bbox,
                           nBlocks = x$nBlocks, nTiles = length(x$tileTable))
    cat("Overall number of events:", x$nEvents,
        if (x$nTypes==1) "(single type)" else paste0("(",x$nTypes," types)"),
        "\n")
    
    cat("\nSummary of event marks and number of potential sources:\n")
    print(summary(cbind(x$eventMarks, "|.sources|"=x$nSources)), ...)

    invisible(x)
}

as.stepfun.epidataCS <- function (x, ...)
{
    eventTimes <- x$events$time
    removalTimes <- eventTimes + x$events$eps.t
    tps <- sort(unique(c(eventTimes, removalTimes[is.finite(removalTimes)])))
    nInfectious <- sapply(tps, function(t) sum(eventTimes <= t & removalTimes > t))
    stepfun(tps, c(0,nInfectious), right = TRUE) # no ties, 'tps' is unique
}



###################################################
### Distances from potential (eps.s, eps.t) sources
###################################################

getSourceDists <- function (object, dimension = c("space", "time"))
{
    dimension <- match.arg(dimension)
    
    ## extract required info from "epidataCS"-object
    distmat <- as.matrix(dist(
        if (dimension == "space") {
            coordinates(object$events)
        } else object$events$time
        ))
    .sources <- object$events$.sources

    ## number of sources
    nsources <- sapply(.sources, length)
    hasSources <- nsources > 0
    cnsources <- c(0, cumsum(nsources))

    ## generate vector of distances of events to their potential sources
    sourcedists <- numeric(sum(nsources))
    for (i in which(hasSources)) {
        .sourcesi <- .sources[[i]]
        .sourcedists <- distmat[i, .sourcesi]
        .idx <- cnsources[i] + seq_len(nsources[i])
        sourcedists[.idx] <- .sourcedists
        names(sourcedists)[.idx] <- paste(i, .sourcesi, sep="<-")
    }

    ## Done
    sourcedists
}
