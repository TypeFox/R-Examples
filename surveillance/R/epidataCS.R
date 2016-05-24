################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Data structure for CONTINUOUS SPATIO-temporal infectious disease case data
### and a spatio-temporal grid of endemic covariates
###
### Copyright (C) 2009-2015 Sebastian Meyer
### $Revision: 1163 $
### $Date: 2015-01-12 17:42:15 +0100 (Mon, 12. Jan 2015) $
################################################################################


######################################################################
# MAIN GENERATOR FUNCTION FOR epidataCS OBJECTS
# PARAMS:
# events: SpatialPointsDataFrame of cases with obligatory columns
#   time: time point of event
#   tile: reference to spatial unit (tile) in stgrid, where the event is located
#   type: optional type of event (-> marked twinstim). will be converted to a factor variable.
#   eps.t: maximal temporal influence radius (e.g. length of infectious period, time to culling, etc.), may be Inf
#   eps.s: maximal spatial influence radius (e.g. 100 [km]), may be Inf
#   The remaining columns are further marks of the event, e.g. sex, age of infected person (-> epidemic covariates)
#   The column names ".obsInfLength", ".bdist", ".influenceRegion", and ".sources" are reserved.
#   ".obsInfLength": observed length of the infectious period (being part [0,T])
#   ".bdist": minimal distance of the event locations to the boundary
#   ".influenceRegion": object of class "owin", the intersection of W with b(s,eps.s), with origin at s
#   ".sources": potential sources of infection
# stgrid: data.frame with obligatory columns
#   tile: ID of spatial unit (e.g. id of municipality)
#   start, stop: temporal interval
#   area: area of the spatial unit (tile)
#   The remaining columns are endemic covariates.
#   The column name "BLOCK" is reserved (indexing the time intervals of stgrid).
# W: SpatialPolygons. Observation region. Must have same proj4string as events.
# qmatrix: square indicator matrix (0/1 or TRUE/FALSE) for possible transmission between the event types. will be internally converted to logical. Defaults to an independent spread of the event types.
# nCircle2Poly: accuracy (number of edges) of the polygonal approximation of a circle
# T: end of observation period (=last stop time). Must be specified if only the
#    start but not the stop times are supplied in stgrid (-> auto-generation of stop-times).
# clipper: engine to use for computing polygon intersections.
######################################################################

obligColsNames_events <- c("time", "tile", "type", "eps.t", "eps.s")
obligColsNames_stgrid <- c("start", "stop", "tile", "area")
reservedColsNames_events <- c(".obsInfLength", ".sources", ".bdist",
                              ".influenceRegion", "BLOCK", "start")
reservedColsNames_stgrid <- c("BLOCK")

as.epidataCS <- function (events, stgrid, W, qmatrix = diag(nTypes),
                          nCircle2Poly = 32, T = NULL,
                          clipper = c("polyclip", "rgeos"),
                          verbose = interactive())
{
    clipper <- match.arg(clipper)
    
    # Check and SORT events
    if (verbose) cat("\nChecking 'events':\n")
    events <- check_events(events, verbose = verbose)

    # Check and SORT stgrid
    if (verbose) cat("Checking 'stgrid':\n")
    tiles <- NULL                       # FIXME: add argument to as.epidataCS
    stgrid <- if (missing(stgrid) && inherits(tiles, "SpatialPolygons")) {
        if (verbose) cat("\t(missing, using time-constant 'tiles' grid)\n")
        stgrid_template <- data.frame(
            start = min(events$time),
            stop = max(events$time),
            tile = row.names(tiles),
            area = areaSpatialPolygons(tiles, byid = TRUE),
            check.rows = FALSE, check.names = FALSE)
        check_stgrid(stgrid_template, verbose = FALSE)
    } else {
        check_stgrid(stgrid, T, verbose = verbose)
    }
    
    # Check class of W and consistency of area
    if (verbose) cat("Checking 'W' ...\n")
    W <- check_W(W, area.other =
                 sum(stgrid[["area"]][seq_len(nlevels(stgrid$tile))]),
                 other = "stgrid")
    stopifnot(identicalCRS(W, events))
    
    # Check qmatrix
    if (verbose) cat("Checking 'qmatrix' ...\n")
    typeNames <- levels(events$type)
    nTypes <- length(typeNames)     # default value of qmatrix depends on nTypes
    qmatrix <- checkQ(qmatrix, typeNames)

    # Check nCircle2Poly
    stopifnot(isScalar(nCircle2Poly))
    nCircle2Poly <- as.integer(nCircle2Poly)
    
    # Small helper function converting event index to (time, tile, type) string
    eventidx2string <- function (eventIdx) {
        with(events@data,
             paste(c("time", "tile", "type"), "=",
                   c(time[eventIdx], dQuote(tile[eventIdx]),
                     dQuote(type[eventIdx])),
                   collapse = ", "))
    }

    # Check that all events are part of W
    if (verbose) cat("Checking if all events are part of 'W' ...\n")
    WIdxOfEvents <- over(events, W)
    if (eventNotInWidx <- match(NA, WIdxOfEvents, nomatch = 0L)) {
        stop("the event at (", eventidx2string(eventNotInWidx), ") is not ",
             "inside 'W'")
    }

    # Some basic quantities
    eventCoords <- coordinates(events)
    nEvents <- nrow(eventCoords)
    timeRange <- with(stgrid, c(start[1], stop[length(stop)]))

    # Are event times covered by stgrid?
    if (verbose) cat("Checking if all events are covered by 'stgrid' ...\n")
    ## FIXME: what about pre-history events? don't need stgrid-data for them
    if (events$time[1L] <= timeRange[1L] || events$time[nEvents] > timeRange[2L]) {
        stop("event times are not covered by 'stgrid': must be in (",
             timeRange[1L],",",timeRange[2L],"]")
    }

    # Are all events$tile references really part of the stgrid?
    .events.tile <- factor(events$tile, levels = levels(stgrid$tile))
    if (missingSCellIdx <- match(NA, .events.tile, nomatch = 0L)) {
        stop("the 'events$tile' entry \"", events$tile[missingSCellIdx], "\"",
             " is not a valid level of 'stgrid$tile'")
    }
    events$tile <- .events.tile

    # Calculate time point of removal, when event is definitely no longer infective
    removalTimes <- events$time + events$eps.t

    # Calculate distance matrix of events
    if (verbose) cat("Calculating Euclidean distance matrix of events ...\n")
    eventDists <- as.matrix(dist(eventCoords, method = "euclidean"))
    #diag(eventDists) <- Inf   # infinite distance to oneself (no self-infection), not needed

    # Map events to corresponding grid cells
    # Also precalculate possible origins of events (other infected individuals)
    if (verbose) cat("Mapping events to 'stgrid' cells and",
                     "determining potential event sources ...\n")
    withPB <- verbose && interactive()
    gridcellsOfEvents <- integer(nEvents)
    eventSources <- vector(nEvents, mode = "list")
    if (withPB) pb <- txtProgressBar(min=0, max=nEvents, initial=0, style=3)
    for (i in seq_len(nEvents)) {
        idx <- gridcellOfEvent(events$time[i], events$tile[i], stgrid)
        if (is.na(idx)) {
            stop("could not find information for time point ", events$time[i],
                 " and tile \"", events$tile[i], "\" in 'stgrid'")
        }
        gridcellsOfEvents[i] <- idx
        eventSources[[i]] <- determineSources(
            i, events$time, removalTimes, eventDists[i,], events$eps.s, events$type, qmatrix
        )
        if (withPB) setTxtProgressBar(pb, i)
    }
    if (withPB) close(pb)

    # Attach endemic covariates from stgrid to events
    if (verbose) cat("Attaching endemic covariates from 'stgrid' to 'events' ...\n")
    stgridIgnoreCols <- match(setdiff(obligColsNames_stgrid, "start"), names(stgrid))
    copyCols <- setdiff(seq_along(stgrid), stgridIgnoreCols)
    reservedColsIdx <- na.omit(match(names(stgrid)[copyCols], names(events@data),
                                     nomatch=NA_integer_))
    if (length(reservedColsIdx) > 0L) {
        warning("in 'events@data', the existing columns with names of endemic ",
                "covariates from 'stgrid' (",
                paste0("'", names(events@data)[reservedColsIdx], "'", collapse=", "),
                ") have been replaced")
        events@data <- events@data[-reservedColsIdx]
    }
    events@data <- cbind(events@data, stgrid[gridcellsOfEvents, copyCols])

    # Calculate observed infection length = min(T-time, eps.t) for use in log-likelihood
    events$.obsInfLength <- with(events@data, pmin(timeRange[2]-time, eps.t))

    # Attach possible eventSources (infective individuals) to events
    events$.sources <- eventSources

    # Calculate minimal distance of event locations from the polygonal boundary
    if (verbose) cat("Calculating the events' distances to the boundary ...\n")
    Wowin <- as(W, "owin")              # imported from polyCub
    events$.bdist <- bdist(eventCoords, Wowin)

    # Construct spatial influence regions around events
    if (verbose) cat("Constructing spatial influence regions around events ...\n")
    events$.influenceRegion <- if (clipper == "polyclip") {
        .influenceRegions(events, Wowin, nCircle2Poly, clipper=clipper)
    } else .influenceRegions(events, W, nCircle2Poly, clipper=clipper)

    # Return components in a list of class "epidataCS"
    res <- list(events = events, stgrid = stgrid, W = W, qmatrix = qmatrix)
    class(res) <- "epidataCS"
    if (verbose) cat("Done.\n\n")
    return(res)
}






######################################################################
# HELPER FUNCTIONS FOR as.epidataCS
######################################################################


### CHECK FUNCTION FOR events ARGUMENT IN as.epidataCS

check_events <- function (events, dropTypes = TRUE, verbose = TRUE)
{
    # Check class and spatial dimensions
    stopifnot(inherits(events, "SpatialPointsDataFrame"))
    if (ncol(events@coords) != 2L) {
        stop("only two spatial dimensions are supported")
    }

    # check suitability of Euclidean geometry
    if (identical(FALSE, is.projected(events))) { # is.projected may return NA
        warning("\"epidataCS\" expects planar coordinates; ",
                "see 'spTransform' in package \"rgdal\"")
    }

    # Check existence of type column
    if (verbose) cat("\tChecking 'type' column ... ")
    events$type <- if ("type" %in% names(events)) {
                       if (dropTypes) factor(events$type) else as.factor(events$type)
                   } else {
                       if (verbose) cat("Setting 'type' to 1 for all events.")
                       factor(rep.int(1L,nrow(events@coords)))
                     }
    if (verbose) cat("\n")

    # Check obligatory columns
    obligColsIdx <- match(obligColsNames_events, names(events), nomatch = NA_integer_)
    if (any(obligColsMissing <- is.na(obligColsIdx))) {
        stop("missing obligatory columns in 'events@data': ",
            paste(obligColsNames_events[obligColsMissing], collapse = ", "))
    }

    # Check other columns on reserved names
    reservedColsIdx <- na.omit(match(reservedColsNames_events, names(events),
                                     nomatch=NA_integer_))
    if (length(reservedColsIdx) > 0L) {
        warning("in 'events@data', the existing columns with reserved names (",
                paste0("'", names(events)[reservedColsIdx], "'", collapse=", "),
                ") have been replaced")
        events@data <- events@data[-reservedColsIdx]
    }    

    # Check that influence radii are numeric and positive (also not NA)
    if (verbose) cat("\tChecking 'eps.t' and 'eps.s' columns ...\n")
    with(events@data, stopifnot(is.numeric(eps.t), eps.t > 0,
                                is.numeric(eps.s), eps.s > 0))

    # Transform time into a numeric variable
    if (verbose) cat("\tConverting event time into a numeric variable ...\n")
    events$time <- as.numeric(events$time)
    stopifnot(!is.na(events$time))

    # Check event times for ties
    if (verbose) cat("\tChecking event times for ties ...\n")
    timeIsDuplicated <- duplicated(events$time)
    if (any(timeIsDuplicated)) {
        duplicatedTimes <- unique(events$time[timeIsDuplicated])
        warning("detected non-unique event times: ",
                "concurrent events at time ",
                if (length(duplicatedTimes) == 1L) "point " else "points\n",
                paste(duplicatedTimes, collapse = ", "))
    }

    # Sort events chronologically
    if (verbose) cat("\tSorting events ...\n")
    events <- events[order(events$time),]
    
    # First obligatory columns then remainders (epidemic covariates)
    obligColsIdx <- match(obligColsNames_events, names(events@data))
    covarColsIdx <- setdiff(seq_along(events@data), obligColsIdx)
    events@data <- events@data[c(obligColsIdx, covarColsIdx)]
    events@coords.nrs <- numeric(0L)  # forget index of coordinate columns

    # Done.
    return(events)
}



### CHECK FUNCTION FOR stgrid ARGUMENT IN as.epidataCS

check_stgrid <- function (stgrid, T, verbose = TRUE)
{
    # Check class
    stopifnot(inherits(stgrid, "data.frame"))

    # Check obligatory columns
    autostop <- FALSE
    if (is.null(stgrid[["stop"]])) {
        if (is.null(T)) stop("'T' must be specified for auto-generation ",
                             "of 'stop' column in 'stgrid'")
        stopifnot(isScalar(T))
        autostop <- TRUE
        stgrid$stop <- NA_real_
    }
    obligColsIdx <- match(obligColsNames_stgrid, names(stgrid), nomatch = NA_integer_)
    if (any(obligColsMissing <- is.na(obligColsIdx))) {
        stop("missing obligatory columns in 'stgrid': ",
            paste(obligColsNames_stgrid[obligColsMissing], collapse = ", "))
    }

    # Check other columns on reserved names
    reservedColsIdx <- na.omit(match(reservedColsNames_stgrid, names(stgrid),
                                     nomatch=NA_integer_))
    if (length(reservedColsIdx) > 0L) {
        warning("in 'stgrid', the existing columns with reserved names (",
                paste0("'", names(stgrid)[reservedColsIdx], "'", collapse=", "),
                ") have been replaced")
        stgrid <- stgrid[-reservedColsIdx]
    }

    # Transform tile into a factor variable
    # (also removing unused levels if it was a factor)
    if (verbose) cat("\tConverting 'tile' into a factor variable ...\n")
    stgrid$tile <- factor(stgrid$tile)

    # Transform start times and area into numeric variables
    stgrid$start <- as.numeric(stgrid$start)
    stgrid$area <- as.numeric(stgrid$area)        

    # Check stop times
    stgrid$stop <- if (autostop) {
        # auto-generate stop times from start times and T
        if (verbose) cat("\tAuto-generating 'stop' column ...\n")
        starts <- sort(unique(stgrid$start))
        if (T <= starts[length(starts)]) {
            stop("'T' must be larger than the last 'start' time in 'stgrid'")
        }
        stops <- c(starts[-1], T)
        stops[match(stgrid$start, starts)]
    } else {
        as.numeric(stgrid$stop)
    }

    # chronological data.frame of unique periods
    histIntervals <- unique(stgrid[c("start", "stop")])
    histIntervals <- histIntervals[order(histIntervals[,1L]),]
    nBlocks <- nrow(histIntervals)

    if (!autostop) {
        # Check start/stop consistency
        if (verbose) cat("\tChecking start/stop consisteny ...\n")
        if (any(histIntervals[,2L] <= histIntervals[,1L])) {
            stop("stop times must be greater than start times")
        }
        startStopCheck <- histIntervals[-1L,1L] != histIntervals[-nBlocks,2L]
        if (startStopCheckIdx <- match(TRUE, startStopCheck, nomatch = 0)) {
            stop("inconsistent start/stop times: time intervals not consecutive ",
                 "at stop time ", histIntervals[startStopCheckIdx,2L])
        }
    }

    # Add BLOCK id
    stgrid$BLOCK <- match(stgrid$start, histIntervals[,1L])

    # Check that we have a full BLOCK x tile grid
    if (verbose) cat("\tChecking if the grid is complete ...\n")
    blocksizes <- table(stgrid$BLOCK)
    tiletable <- table(stgrid$tile)
    if (length(unique(blocksizes)) > 1L || length(unique(tiletable)) > 1L) {
        stop("'stgrid' is not a full grid")
    }

    # First column BLOCK, then obligCols, then remainders (endemic covariates)
    if (verbose) cat("\tSorting the grid by time and tile ...\n")
    BLOCKcolIdx <- match("BLOCK", names(stgrid))
    obligColsIdx <- match(obligColsNames_stgrid, names(stgrid))
    covarColsIdx <- setdiff(seq_along(stgrid), c(BLOCKcolIdx, obligColsIdx))
    stgrid <- stgrid[c(BLOCKcolIdx, obligColsIdx, covarColsIdx)]

    # Sort by BLOCK and tile
    stgrid <- stgrid[order(stgrid$BLOCK, stgrid$tile),]

#     # Get row indexes of the blocks' first/last rows
#     beginBlock <- match(seq_len(nBlocks), stgrid[["BLOCK"]])
#     endBlock <- c(beginBlock[-1L]-1L, nrow(stgrid))

    # Done.
    return(stgrid)
}



### CHECK FUNCTION FOR W ARGUMENT IN as.epidataCS

check_W <- function (W, area.other = NULL, other, tolerance = 0.001)
{
    W <- as(W, "SpatialPolygons") # i.e. drop data if a SpatialPolygonsDataFrame
    
    if (!is.null(area.other) && area.other > 0) {
        check_W_area(W, area.other, other, tolerance)
    }
    
    return(W)
}

check_W_area <- function (W, area.other, other, tolerance = 0.001)
{
    area.W <- areaSpatialPolygons(W)
    if (!isTRUE(all.equal.numeric(area.other, area.W, tolerance = tolerance,
                                  check.attributes = FALSE)))
        warning("area of 'W' (", area.W, ") differs from ",
                "total tile area in '", other, "' (", area.other, ")")
}



### CHECK FUNCTION FOR tiles ARGUMENT IN as.epidataCS

check_tiles <- function (tiles, levels,
                         events = NULL, areas.stgrid = NULL, W = NULL,
                         keep.data = FALSE, tolerance = 0.05)
{
    stopifnot(inherits(tiles, "SpatialPolygons"),
              is.vector(levels, mode="character"))
    tileIDs <- row.names(tiles)

    ## check completeness of tiles
    if (any(missingtiles <- !levels %in% tileIDs))
        stop(sum(missingtiles), " regions are missing in 'tiles', ",
             "check 'row.names(tiles)'")

    ## re-order: first 'levels', then any extra tiles
    tiles <- tiles[c(levels, setdiff(tileIDs, levels)),]

    ## drop data (also for suitable over-method in check_tiles_events)
    .tiles <- as(tiles, "SpatialPolygons")
    
    ## check tile specification of events and identical projection
    if (!is.null(events)) {
        check_tiles_events(.tiles, events)
    }

    ## check areas
    areas.tiles <- areaSpatialPolygons(tiles[levels,], byid = TRUE)
    if (!is.null(areas.stgrid)) {
        check_tiles_areas(areas.tiles, areas.stgrid, tolerance=tolerance)
    }
    if (!is.null(W)) {
        stopifnot(identicalCRS(tiles, W))
        check_W_area(W, area.other=sum(areas.tiles), other="tiles",
                     tolerance=tolerance)
    }

    ## done
    if (keep.data) tiles else .tiles
}

check_tiles_events <- function (tiles, events)
{
    tiles <- as(tiles, "SpatialPolygons") # remove potential data for over()
    stopifnot(inherits(events, "SpatialPointsDataFrame"),
              identicalCRS(tiles, events))
    tileIDs <- row.names(tiles)
    eventIDs <- row.names(events)
    
    ## get polygon ID's of events (via overlay)
    eventtiles <- tileIDs[over(events, tiles)]
    
    if (length(which_not_in_tiles <- which(is.na(eventtiles))))
        warning("some of 'events' are not within 'tiles': ",
                paste0("\"", eventIDs[which_not_in_tiles], "\"", collapse=", "))

    if (!is.null(events@data[["tile"]])) {
        which_disagree <- setdiff(
            which(eventtiles != as.character(events$tile)),
            which_not_in_tiles)
        if (length(which_disagree))
            message("'over(events, tiles)' disagrees with 'events$tile': ",
                    paste0("\"", eventIDs[which_disagree], "\"", collapse=", "))
    }
    invisible()
}

check_tiles_areas <- function (areas.tiles, areas.stgrid, tolerance = 0.05)
{
    areas_all_equal <- all.equal.numeric(areas.stgrid, areas.tiles,
                                         tolerance = tolerance,
                                         check.attributes = FALSE)
    if (!isTRUE(areas_all_equal))
        warning("tile areas in 'stgrid' differ from areas of 'tiles': ",
                areas_all_equal)
}


### CONSTRUCT SPATIAL INFLUENCE REGIONS AROUND EVENTS

# An influenceRegion is an object of class "owin" with origin
# at the event (over which we have to integrate by a cubature rule)
# An attribute "area" gives the area of the influenceRegion.
# If it is actually a circular influence region, then there is an attribute
# "radius" denoting the radius of the influence region.
# Argument 'W' can be of class "owin" (preferred) or "SpatialPolygons"
# (especially for clipper="rgeos")
.influenceRegions <- function (events, W, npoly, maxExtent = NULL,
                               clipper = "polyclip")
{
    Wowin <- as(W, "owin")
    if (is.null(maxExtent)) maxExtent <- diameter.owin(Wowin)
    doIntersection <- switch(
        clipper,  # which package to use for polygon intersection
        "polyclip" = function (center, eps)
            intersectPolyCircle.owin(Wowin, center, eps, npoly),
        "rgeos" = function (center, eps) as(
            intersectPolyCircle.SpatialPolygons(
                as(W, "SpatialPolygons"), center, eps, npoly),
            "owin"),
        stop("unsupported polygon clipping engine: '", clipper, "'")
        )
    
    eventCoords <- coordinates(events)
    res <- mapply(
        function (x, y, eps, bdist) {
            center <- c(x,y)
            ## if eps is very large, the influence region is the whole region of W
            iR <- shift.owin(
                if (eps > maxExtent) Wowin else doIntersection(center, eps),
                -center)
            ## if iR is actually a circle of radius eps, attach eps as attribute
            attr(iR, "area") <- if (eps <= bdist) {
                attr(iR, "radius") <- eps
                pi * eps^2
            } else area.owin(iR)
            iR
        },
        eventCoords[,1], eventCoords[,2], events$eps.s, events$.bdist,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)
    
    attr(res, "nCircle2Poly") <- npoly
    attr(res, "clipper") <- clipper
    res
}
