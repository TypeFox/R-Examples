################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Data structure "epidata" representing the SIR event history of a fixed
### geo-referenced population (e.g., farms, households) for twinSIR() analysis
###
### Copyright (C) 2008-2010, 2012, 2014-2015 Sebastian Meyer
### $Revision: 1518 $
### $Date: 2015-11-13 13:22:31 +0100 (Fre, 13. Nov 2015) $
################################################################################

## CAVE:
## - we assume fixed coordinates (this is important since time-varying
##   coordinates would result in more sophisticated and time consuming
##   calculations of distance matrices) !
## - in the first block (start = t0) all id's must be present (for coordinates)
## - those id's with atRiskY(t0) = 0 are taken as initially infectious
## - SIS epidemics are possible, but must be given as SIRS with pseudo R-events,
##   i.e. individuals will be removed and become susceptible directly afterwards


################################################################################
## Convert a simple data.frame with one row per individual and with columns for
## the times of becoming exposed/infectious/removed
## to the long "epidata" event history start/stop format.
## tE.col and tR.col can be missing corresponding to SIR, SEI, or SI data.
## NA's in time variables mean that the respective event has not yet occurred.
## Time-varying covariates are not supported by this converter.
################################################################################

as.epidata.data.frame <- function (data, t0, tE.col, tI.col, tR.col,
                                   id.col, coords.cols, f = list(), w = list(),
                                   D = dist, keep.cols = TRUE, ...)
{
    if (missing(t0)) {
        return(NextMethod("as.epidata"))  # as.epidata.default
    }

    ## drop individuals that have already been removed prior to t0
    ## since they would otherwise be considered as initially infective
    ## (atRiskY = 0 in first time block) and never be removed
    if (!missing(tR.col)) {
        alreadyRemoved <- !is.na(data[[tR.col]]) & data[[tR.col]] <= t0
        if (any(alreadyRemoved)) {
            data <- data[!alreadyRemoved,]
            message("Note: dropped rows with tR <= t0 (",
                    paste0(which(alreadyRemoved), collapse = ", "), ")")
        }
    }

    ## parse id column
    id <- factor(data[[id.col]]) # removes unused levels
    stopifnot(!anyDuplicated(id), !is.na(id))
    N <- nlevels(id) # = nrow(data)

    ## make time relative to t0
    subtract_t0 <- function (x) as.numeric(x - t0)
    tI <- subtract_t0(data[[tI.col]])
    tE <- if (missing(tE.col)) tI else subtract_t0(data[[tE.col]])
    tR <- if (missing(tR.col)) rep.int(NA_real_, N) else subtract_t0(data[[tR.col]])

    ## check E-I-R order
    if (any((is.na(tE) & !(is.na(tI) & is.na(tR))) | (is.na(tI) & !is.na(tR)))) {
        stop("events cannot be skipped (NA in E/I => NA in I/R)")
    }
    if (any(.wrongsequence <- (tE > tI | tI >= tR) %in% TRUE)) {  # TRUE | NA = TRUE
        stop("E-I-R events are in wrong order for the following id's: ",
             paste0(id[.wrongsequence], collapse = ", "))
    }

    ## vector of stop times
    stopTimes <- c(tE, tI, tR)
    stopTimes <- stopTimes[!is.na(stopTimes) & stopTimes > 0]
    stopTimes <- sort.int(unique.default(stopTimes), decreasing = FALSE)
    nBlocks <- length(stopTimes)
    if (nBlocks == 0L) {
        stop("nothing happens after 't0'")
    }
    
    ## initialize event history
    evHist <- data.frame(
        id = rep.int(id, nBlocks),
        start = rep.int(c(0,stopTimes[-nBlocks]), rep.int(N, nBlocks)),
        stop = rep.int(stopTimes, rep.int(N, nBlocks)),
        atRiskY = NA, event = 0, Revent = 0, # adjusted in the loop below
        row.names = NULL, check.rows = FALSE, check.names = FALSE)
    
    ## indexes of the last rows of the time blocks
    blockbase <- c(0, seq_len(nBlocks) * N)

    ## which individuals are at risk in the first (next) block
    Y <- is.na(tE) | tE > 0
    
    ## Loop over the blocks/stop times to adjust atRiskY, event and Revent
    for (i in seq_len(nBlocks)) {
        ct <- stopTimes[i]

        ## set individual at-risk indicators for the current time block
        evHist$atRiskY[blockbase[i] + seq_len(N)] <- Y
        ## individuals who become exposed at the current stop time
        ## will no longer be at risk in the next block
        Y[which(tE == ct)] <- FALSE
        
        ## process events at this stop time
        evHist$event[blockbase[i] + which(tI == ct)] <- 1
        evHist$Revent[blockbase[i] + which(tR == ct)] <- 1
    }

    ## add additional time-constant covariates
    extraVarNames <- coords.cols  # may be NULL
    if (isTRUE(keep.cols)) {
        extraVarNames <- c(extraVarNames, setdiff(names(data), id.col))
    } else if (length(keep.cols) > 0L && !identical(FALSE, keep.cols)) {
        extraVarNames <- c(extraVarNames, names(data[keep.cols]))
    }
    extraVarNames <- unique.default(extraVarNames)
    if (length(extraVarNames) > 0L) {
        evHist <- data.frame(
            evHist,
            data[rep.int(seq_len(N), nBlocks), extraVarNames, drop=FALSE],
            row.names = NULL, check.names = TRUE, stringsAsFactors = TRUE)
    }

    ## Now we can pass the generated event history to the default method
    ## for the usual consistency checks and the pre-calculation of f covariates
    as.epidata.default(
        data = evHist,
        id.col = "id", start.col = "start", stop.col = "stop",
        atRiskY.col = "atRiskY", event.col = "event", Revent.col = "Revent",
        coords.cols = coords.cols, f = f, w = w, D = D)
}


################################################################################
# DEFAULT CONVERTER, which requires a start/stop event history data.frame
# It performs consistency checks, and pre-calculates the distance-based
# epidemic covariates from f.
################################################################################

as.epidata.default <- function(data, id.col, start.col, stop.col, atRiskY.col,
    event.col, Revent.col, coords.cols, f = list(), w = list(), D = dist, ...)
{
    cl <- match.call()
    
    # If necessary, convert 'data' into a data.frame (also converting
    # column names to syntactically correct names for use in formulae)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
        
    # Use column numbers as indices and check them
    colargs <- c("id.col", "start.col", "stop.col", "atRiskY.col",
                 "event.col", "Revent.col", "coords.cols")
    colidxs <- structure(as.list(numeric(length(colargs))), names = colargs)
    for (colarg in colargs) {
        colidx <- get(colarg, inherits = FALSE)
        if (colarg != "coords.cols" && length(colidx) != 1L) {
            stop("the column specifier '", colarg, "' must be of length 1")
        }
        if (is.character(colidx)) {
            colidx <- match(colidx, colnames(data))
            if (any(is.na(colidx))) {
                stop("'", colarg, " = ", deparse(cl[[colarg]]), "': ",
                     "column does not exist in 'data'")
            }
        } else if (is.numeric(colidx) && any(colidx<1L | colidx>ncol(data))) {
            stop("'", colarg, " = ", deparse(cl[[colarg]]), "': ",
                 "column index must be in [1; ", ncol(data), "=ncol(data)]")
        }
        colidxs[[colarg]] <- colidx
    }
    
    # Rename main columns to default column names
    colidxsVec <- unlist(colidxs)
    colnams <- c("id", "start", "stop", "atRiskY", "event", "Revent")
    colnames(data)[colidxsVec[1:6]] <- colnams
    usedReservedName <- any(colnams %in% colnames(data)[-colidxsVec[1:6]])
    
    # REORDER COLUMNS, so that main columns come first (also for make.unique)
    data <- data[c(colidxsVec, setdiff(seq_len(NCOL(data)), colidxsVec))]
    
    # Make columns names unique (necessary if other column with name in colnams)
    if (usedReservedName) {
        colnames(data) <- make.unique(colnames(data))
        message("Some other columns had reserved names and have been renamed")
    }
    
    # Convert id into a factor (also removing unused levels if it was a factor)
    data[["id"]] <- factor(data[["id"]])
    
    # Check atRiskY, event and Revent for values other than 0 and 1
    for (var in c("atRiskY", "event", "Revent")) {
        data[[var]] <- as.numeric(data[[var]])
        if (any(! data[[var]] %in% c(0,1)))
            stop("'", var, "' column may only assume values 0 and 1")
    }
    
    # Check consistency of atRiskY and event (event only if at-risk)
    noRiskButEvent <- data[["atRiskY"]] == 0 & data[["event"]] == 1
    if (noRiskButEventRow <- match(TRUE, noRiskButEvent, nomatch = 0)) {
        stop("inconsistent atRiskY/event indicators in row ",
             noRiskButEventRow, ": event only if at risk")
    }

    # Check event (infection) times for ties
    eventTimes <- data[data[["event"]] == 1, "stop"]
    ReventTimes <- data[data[["Revent"]] == 1, "stop"]
    duplicatedEventTime <- duplicated(c(eventTimes, ReventTimes))
    if (duplicatedEventTimeIdx <- match(TRUE, duplicatedEventTime, nomatch=0)) {
        stop("non-unique event times: concurrent event/Revent at time ",
             c(eventTimes, ReventTimes)[duplicatedEventTimeIdx])
    }

    # Check start/stop consistency and add block id
    histIntervals <- unique(data[c("start", "stop")])
    histIntervals <- histIntervals[order(histIntervals[,1L]),]
    nBlocks <- nrow(histIntervals)
    if (any(histIntervals[,2L] <= histIntervals[,1L])) {
        stop("stop times must be greater than start times")
    }
    startStopCheck <- histIntervals[-1L,1L] != histIntervals[-nBlocks,2L]
    if (startStopCheckIdx <- match(TRUE, startStopCheck, nomatch = 0)) {
        stop("inconsistent start/stop times: time intervals not consecutive ",
             "at stop time ", histIntervals[startStopCheckIdx,2L])
    }
    if ("BLOCK" %in% colnames(data)) {
        warning("column name 'BLOCK' is reserved, ",
                "existing column has been replaced")
    }
    data[["BLOCK"]] <- match(data[["start"]], histIntervals[,1L])
    
    # SORT by block/id and create indexes for block borders
    data <- data[order(data[["BLOCK"]], data[["id"]]),]
    beginBlock <- match(seq_len(nBlocks), data[["BLOCK"]])
    endBlock <- c(beginBlock[-1L]-1L, nrow(data))
    
    # make block column the first column
    BLOCK.col <- match("BLOCK", colnames(data))
    data <- data[c(BLOCK.col, setdiff(seq_along(data), BLOCK.col))]
    coords.cols <- 1L + 6L + seq_along(colidxs[["coords.cols"]])

    # Check consistency of atRiskY and event (not at-risk after event) 
    .checkFunction <- function(eventblock, eventid)
    {
        if (eventblock == nBlocks) return(invisible())
        rowsOfNextBlock <- beginBlock[eventblock+1L]:endBlock[eventblock+1L]
        nextBlockData <- data[rowsOfNextBlock, c("id", "atRiskY")]
        idIdx <- which(nextBlockData[["id"]] == eventid)
        if (length(idIdx) == 1L && nextBlockData[idIdx, "atRiskY"] == 1) {
            stop("inconsistent atRiskY/event indicators for id '", eventid,
                 "': should not be at risk immediately after event")
        }
    }
    eventTable <- data[data[["event"]] == 1,]
    for(k in seq_len(nrow(eventTable)))
    {
        .checkFunction(eventTable[k,"BLOCK"], eventTable[k,"id"])
    }

    # Set attributes
    attr(data, "eventTimes") <- sort(eventTimes)
    attr(data, "timeRange") <- c(histIntervals[1L,1L],histIntervals[nBlocks,2L])
    attr(data, "coords.cols") <- coords.cols
    # <- must include this info because externally of this function
    #    we don't know how many coords.cols (dimensions) we have
    attr(data, "f") <- list()  # initialize
    attr(data, "w") <- list()  # initialize
    class(data) <- c("epidata", "data.frame")

    # Compute epidemic variables
    update.epidata(data, f = f, w = w, D = D)
}


update.epidata <- function (object, f = list(), w = list(), D = dist, ...)
{
    oldclass <- class(object)
    class(object) <- "data.frame" # avoid use of [.epidata

    ## block indexes and first block
    beginBlock <- which(!duplicated(object[["BLOCK"]],
                                    nmax = object[["BLOCK"]][nrow(object)]))
    endBlock <- c(beginBlock[-1L]-1L, nrow(object))
    firstDataBlock <- object[seq_len(endBlock[1L]), ]

    ## check f and calculate distance matrix
    if (length(f) > 0L) {
        if (!is.list(f) || is.null(names(f)) || any(!sapply(f, is.function))) {
            stop("'f' must be a named list of functions")
        }
        lapply(X = f, FUN = function (B) {
            if (!isTRUE(all.equal(c(5L,2L), dim(B(matrix(0, 5, 2))))))
                stop("'f'unctions must retain the dimensions of their input")
        })
        if (any(names(f) %in% names(object))) {
            warning("'f' components replace existing columns of the same name")
        }

        ## reset / initialize columns for distance-based epidemic weights
        object[names(f)] <- 0
        ## keep functions as attribute
        attr(object, "f")[names(f)] <- f
        
        ## check / compute distance matrix
        distmat <- if (is.function(D)) {
            if (length(coords.cols <- attr(object, "coords.cols")) == 0L) {
                stop("need coordinates to calculate the distance matrix")
            }
            coords <- as.matrix(firstDataBlock[coords.cols],
                                rownames.force = FALSE)
            rownames(coords) <- as.character(firstDataBlock[["id"]])
            as.matrix(D(coords))
        } else { # a numeric matrix (or "Matrix")
            if (length(dn <- dimnames(D)) != 2L) {
                stop("if not a function, 'D' must be a matrix-like object")
            }
            if (!all(firstDataBlock[["id"]] %in% dn[[1L]],
                     firstDataBlock[["id"]] %in% dn[[2L]])) {
                stop("'dimnames(D)' must contain the individuals' IDs")
            }
            D
        }
    }

    ## check covariate-based epidemic weights
    if (length(w) > 0L) {
        if (!is.list(w) || is.null(names(w)) || any(!sapply(w, is.function))) {
            stop("'w' must be a named list of functions")
        }
        if (any(names(w) %in% names(object))) {
            warning("'w' components replace existing columns of the same name")
        }
        
        ## reset / initialize columns for covariate-based epidemic weights
        object[names(w)] <- 0
        ## keep functions as attribute
        attr(object, "w")[names(w)] <- w

        ## compute wij matrix for each of w
        wijlist <- compute_wijlist(w = w, data = firstDataBlock)
    }

    ## Compute sum of epidemic covariates over infectious individuals
    if (length(f) + length(w) > 0L) {
        infectiousIDs <- firstDataBlock[firstDataBlock[["atRiskY"]] == 0, "id"]
        ##<- this is a factor variable
        for(i in seq_along(beginBlock)) {
            blockidx <- beginBlock[i]:endBlock[i]
            blockdata <- object[blockidx,]
            blockIDs <- blockdata[["id"]]
            if (length(infectiousIDs) > 0L) {
                if (length(f) > 0L) {
                    u <- distmat[as.character(blockIDs),
                                 as.character(infectiousIDs),
                                 drop = FALSE] # index by factor levels
                    object[blockidx,names(f)] <- vapply(
                        X = f, FUN = function (B) rowSums(B(u)),
                        FUN.VALUE = numeric(length(blockIDs)),
                        USE.NAMES = FALSE)
                }
                if (length(w) > 0L) {
                    object[blockidx,names(w)] <- vapply(
                        X = wijlist, FUN = function (wij) {
                            ## actually don't have to care about the diagonal:
                            ## i at risk => sum does not include it
                            ## i infectious => atRiskY = 0 (ignored in twinSIR)
                            rowSums(wij[as.character(blockIDs),
                                        as.character(infectiousIDs),
                                        drop = FALSE]) # index by factor levels
                        }, FUN.VALUE = numeric(length(blockIDs)),
                        USE.NAMES = FALSE)
                }
            }
            ## update the set of infectious individuals for the next block
            recoveredID <- blockIDs[blockdata[["Revent"]] == 1]
            infectedID <- blockIDs[blockdata[["event"]] == 1]
            if (length(recoveredID) > 0L) {
                infectiousIDs <- infectiousIDs[infectiousIDs != recoveredID]
            } else if (length(infectedID) > 0L) {
                infectiousIDs[length(infectiousIDs)+1L] <- infectedID
            }
        }
    }
    
    ## restore "epidata" class
    class(object) <- oldclass
    return(object)
}

compute_wijlist <- function (w, data)
{
    ## for each function in 'w', determine the variable on which it acts;
    ## this is derived from the name of the first formal argument, which
    ## must be of the form "varname.i"
    wvars <- vapply(X = w, FUN = function (wFUN) {
        varname.i <- names(formals(wFUN))[[1L]]
        substr(varname.i, 1, nchar(varname.i)-2L)
    }, FUN.VALUE = "", USE.NAMES = TRUE)
    
    if (any(wvarNotFound <- !wvars %in% names(data))) {
        stop("'w' function refers to unknown variables: ",
             paste0(names(w)[wvarNotFound], collapse=", "))
    }
   
    ## compute weight matrices w_ij for each of w
    mapply(
        FUN = function (wFUN, wVAR, ids) {
            wij <- outer(X = wVAR, Y = wVAR, FUN = wFUN)
            dimnames(wij) <- list(ids, ids)
            wij
        },
        wFUN = w, wVAR = data[wvars],
        MoreArgs = list(ids = as.character(data[["id"]])),
        SIMPLIFY = FALSE, USE.NAMES = TRUE
    )
}


################################################################################
# EXTRACTION OPERATOR FOR 'EPIDATA' OBJECTS
# Indexing with "[" would be possible (inheriting from data.frame).
# But using any column index would remove attributes (row indexes would not).
# Thus, we define an own method to retain and adjust the attributes when 
# selecting a subset of blocks of the 'epidata'.
# Selecting a subset of columns will remove class "epidata" (resulting in a
# simple data.frame)
################################################################################

"[.epidata" <- function(x, i, j, drop)
{
    # use data.frame method first
    xx <- NextMethod("[")
    # then return its result as pure data.frame or assure valid 'epidata'
    
    # if a subset of columns has been selected and attributes have been removed
    if (NCOL(xx) != ncol(x) || any(names(xx) != names(x))) {
        if (inherits(xx, "data.frame")) { # xx could be a vector
            class(xx) <- "data.frame"  # remove class 'epidata'
        }
        message("Note: converted class \"epidata\" to simple \"", class(xx),
                "\"")
        return(xx)
    }
    # else there is no effective column selection (e.g. j=TRUE)
    
    if (nrow(xx) == 0) {
        message("Note: no rows selected, dropped class \"epidata\"")
        class(xx) <- "data.frame"
        return(xx[TRUE])   # removes attributes
    }
    
    invalidEpidata <- FALSE
    blocksizesx <- table(x[["BLOCK"]])
    blocksizesxx <- table(xx[["BLOCK"]])
    blocksOK <- identical(c(blocksizesxx), c(blocksizesx[names(blocksizesxx)]))
    if (is.numeric(i) && any(diff(na.omit(i)) < 0)) {
        # epidata should remain ordered by time
        warning("dropped class \"epidata\": reordering rows is not permitted")
        invalidEpidata <- TRUE
    } else if (!blocksOK) {
        # blocks should not be cut, epidemic covariates might become invalid
        warning("dropped class \"epidata\": subsetting blocks not allowed")
        invalidEpidata <- TRUE
    } else if (any(diff(as.numeric(names(blocksizesxx))) != 1)) {
        # blocks can only be selected consecutively
        warning("dropped class \"epidata\": ",
                "only consecutive blocks may be selected")
        invalidEpidata <- TRUE
    }
    
    if (invalidEpidata) {
        class(xx) <- "data.frame"
        xx[TRUE] # removes attributes
    } else {
#         # adjust block index so that it starts at 1
#         firstBlockNumber <- as.numeric(names(blocksizesxx)[1])
#         if (firstBlockNumber > 1) {
#             xx[["BLOCK"]] <- xx[["BLOCK"]] - (firstBlockNumber-1)
#         }
        # Restore or adjust attributes
        tmin <- xx[["start"]][1]
        tmax <- xx[["stop"]][nrow(xx)]
        oldEventTimes <- attr(x, "eventTimes")
        attr(xx, "eventTimes") <-
            if (blocksOK) {
                oldEventTimes[oldEventTimes > tmin & oldEventTimes <= tmax]
            } else {
                xx[["stop"]][xx[["event"]] == 1]
            }
        attr(xx, "timeRange") <- c(tmin, tmax)
        attr(xx, "coords.cols") <- attr(x, "coords.cols")
        attr(xx, "f") <- attr(x, "f")
        xx
    }
}


################################################################################
# INSERT BLOCKS FOR EXTRA STOP TIMES IN 'EPIDATA' OBJECTS
################################################################################

intersperse <- function (epidata, stoptimes, verbose = FALSE)
{
    # Check arguments
    if (!inherits(epidata, "epidata")) {
        stop("'epidata' must inherit from class \"epidata\"")
    }
    if (!is.vector(stoptimes, mode = "numeric")) {
        stop("'stoptimes' must be a numeric vector")
    }
    
    # Identify new 'stoptimes'
    sortedEpiStop <- sort(unique(epidata$stop))
    extraStoptimes <- stoptimes[! stoptimes %in% sortedEpiStop]
    
    # Return original 'epidata' if nothing to do
    if (length(extraStoptimes) == 0) {
#         message("nothing done: no new stop times")
        return(epidata)
    }
    
#    # Retain attributes of 'epidata'
#    .attributes <- attributes(epidata)
#    .attributes <- .attributes[match(c("eventTimes", "timeRange",
#        "coords.cols", "f", "config", "call", "terms"), names(.attributes),
#        nomatch = 0)]

    # Check new 'stoptimes'
    timeRange <- attr(epidata, "timeRange")
    inside <- extraStoptimes > timeRange[1] & extraStoptimes < timeRange[2]
    if (any(!inside)) {
        extraStoptimes <- extraStoptimes[inside]
        warning("ignored extra 'stoptimes' outside the observation period")
    }
    
    # Impute blocks for extraStoptimes
    oldclass <- class(epidata)
    class(epidata) <- "data.frame" # Avoid use of [.epidata (not necessary here)
    blocksize <- sum(epidata$BLOCK == 1)
    nInsert <- length(extraStoptimes)
    lastRow <- nrow(epidata)
    epidata <- rbind(epidata,
                     epidata[rep.int(NA_integer_, nInsert * blocksize),],
                     deparse.level = 0) # add NA rows, to be replaced below
    if (verbose) pb <- txtProgressBar(min=0, max=nInsert, initial=0, style=3)
    for(i in seq_len(nInsert)) {
      extraStop <- extraStoptimes[i]
      nextStoptime <- sortedEpiStop[match(TRUE, sortedEpiStop > extraStop)]
      # Find the block (row indexes) into which the extraStop falls
      rowsMatchedBlock <- which(epidata$stop == nextStoptime)
      # Split this block up into 2 parts
      # later part equals original block with start time = extraStop
      newBlock <- epidata[rowsMatchedBlock,]
      newBlock$start <- extraStop
      # earlier part has stop time = extraStop and no events at this time point
      epidata[rowsMatchedBlock, "stop"] <- extraStop
      epidata[rowsMatchedBlock, "event"] <- 0
      epidata[rowsMatchedBlock, "Revent"] <- 0
      # write the new block to epidata (reorder rows later)
      epidata[lastRow + seq_along(rowsMatchedBlock),] <- newBlock
      lastRow <- lastRow + length(rowsMatchedBlock)
      if (verbose) setTxtProgressBar(pb, i)
    }
    if (verbose) close(pb)
    
    # Adjust BLOCK column
    sortedEpiStop <- sort(c(sortedEpiStop, extraStoptimes))
    epidata$BLOCK <- match(epidata$stop, sortedEpiStop)
    
    # Reorder rows by time and id
    epidata <- epidata[order(epidata$BLOCK, epidata$id), ]
    row.names(epidata) <- NULL
    class(epidata) <- oldclass
    
    return(epidata)
}


################################################################################
# SUMMARY FUNCTION FOR EPIDATA OBJECTS
# the epidemic is summarized by the following returned components:
# - type: one of "SIR", "SI", "SIRS", "SIS"
# - size: number of initially susceptible individuals, which became infected
# - initiallyInfected: vector (factor) of initially infected individuals
# - neverInfected: vector (factor) of never (during the observation period)
#                  infected individuals
# - coordinates: matrix with the coordinates of the individuals (rownames=id's)
# - byID: data.frame with time points of events by id (columns time.I, time.R
#         and optionally time.S)
# - counters: data.frame representing the evolution of the epidemic
################################################################################

summary.epidata <- function (object, ...)
{
    class(object) <- "data.frame"  # avoid use of [.epidata (not necessary here)
    
    # extract coordinates and initially infected individuals
    idlevels <- levels(object[["id"]])
    N <- length(idlevels)
    firstDataBlock <- object[object$BLOCK == min(object$BLOCK),]
    coordinates <- as.matrix(firstDataBlock[attr(object, "coords.cols")])
    rownames(coordinates) <- as.character(firstDataBlock[["id"]])
    initiallyInfected <- firstDataBlock$id[firstDataBlock$atRiskY == 0]
    m <- length(initiallyInfected)
    n <- N - m
    
    ### summary 1: event table with columns id, time and type (of event, S/I/R)
    # Extract time points of the S events for each id
    StimesID <- by(object[c("atRiskY", "stop")], object["id"],
                   function(x) {
                       SeventIdx <- which(diff(x[["atRiskY"]]) == 1)
                       x[["stop"]][SeventIdx]
                   }, simplify=TRUE)
    names(StimesID) <- paste0(names(StimesID), ":")
    StimesVec <- c(unlist(StimesID, use.names = TRUE)) # c() if by() returned an array
    .Sids <- sub("(.+):.*", "\\1", names(StimesVec))
    Stimes <- data.frame(id = factor(.Sids, levels = idlevels),
                         stop = StimesVec, type = rep("S", length(StimesVec)),
                         row.names = NULL, check.names = FALSE,
                         stringsAsFactors = FALSE)
    # Extract time points of the I and R events for each id
    Itimes <- object[object$event == 1, c("id", "stop")]
    Itimes[["type"]] <- rep("I", nrow(Itimes))
    Rtimes <- object[object$Revent == 1, c("id", "stop")]
    Rtimes[["type"]] <- rep("R", nrow(Rtimes))
    
    # Combine the three event types into one data.frame
    eventTable <- rbind(Rtimes, Stimes, Itimes)
      # need this order for the counters below in the case of SIS:
      # pseudo-R-event occures infinitesimally before S
    names(eventTable)[2L] <- "time"
    eventTable <- eventTable[order(eventTable[["id"]], eventTable[["time"]]), ]
    eventTable[["type"]] <- factor(eventTable[["type"]], levels=c("S","I","R"))
    rownames(eventTable) <- NULL
    
    ### summary 2: type and size of the epidemic
    resusceptibility <- length(StimesVec) > 0
    epitype <-
        if (resusceptibility) {
            Rtimes_notLast <- Rtimes[-which.max(Rtimes[,2]),]
            onlyPseudoR <- length(setdiff(Rtimes_notLast[,2], Stimes[,2])) == 0
            if (onlyPseudoR) "SIS" else "SIRS"
        } else {
            if (nrow(Rtimes) > 0) "SIR" else "SI"
        }
    isEverInfected <- idlevels %in% initiallyInfected |
        idlevels %in% unique(eventTable$id[eventTable$type == "I"])
    isNeverInfected <- !isEverInfected
    size <- n - sum(isNeverInfected)
#     everInfected <- factor(idlevels[isEverInfected], levels = idlevels)
    neverInfected <- factor(idlevels[isNeverInfected], levels = idlevels)
    
    ### summary 3: eventTable by id in wide form
    byID_everInfected <-
        if (nrow(eventTable) == 0) {
            data.frame(id = factor(character(0), levels = idlevels),
                       time.I = numeric(0), row.names = NULL,
                       check.names = FALSE, stringsAsFactors = FALSE)
        } else if (!resusceptibility) {
            .res <- reshape(eventTable, direction = "wide", timevar = "type",
                           idvar = "id")
            attr(.res, "reshapeWide") <- NULL
            .res
        } else {
            rowsPerId <- table(eventTable[["id"]])
            modulo3 <- rowsPerId %% 3
            rest1 <- modulo3 == 1
            rest12 <- modulo3 >= 1
            missingR <-
                data.frame(id = names(rowsPerId)[rest1],
                           time = rep(NA_real_, sum(rest1)),
                           type = rep("R", sum(rest1)), row.names = NULL,
                           check.names = FALSE, stringsAsFactors = FALSE)
            missingS <- 
                data.frame(id = names(rowsPerId)[rest12],
                           time = rep(NA_real_, sum(rest12)),
                           type = rep("S", sum(rest12)), row.names = NULL,
                           check.names = FALSE, stringsAsFactors = FALSE)
            eventTable3 <- rbind(eventTable, missingR, missingS)
            eventTable3 <- eventTable3[order(eventTable3[["id"]]),]
            .res <- data.frame(
                eventTable3[eventTable3$type == "I", c("id", "time")],
                eventTable3[eventTable3$type == "R", "time", drop = FALSE],
                eventTable3[eventTable3$type == "S", "time", drop = FALSE],
                row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE
            )
            names(.res) <- c("id", paste("time", c("I", "R", "S"), sep="."))
            .res
        }
    byID_neverInfected <- data.frame(id = neverInfected,
        time.I = rep(NA_real_, n-size), time.R = rep(NA_real_, n-size),
        time.S = rep(NA_real_, n-size), row.names = NULL, check.names = FALSE)
    byID_all <- rbind(byID_everInfected,
                      byID_neverInfected[seq_along(byID_everInfected)])
    byID <- byID_all[order(byID_all[["id"]]),]
    rownames(byID) <- NULL
    
    ### summary 4: upgrade eventTable with
    ###            evolution of numbers of susceptibles, infectious and removed
    counters <- eventTable[order(eventTable[["time"]]),c("time", "type", "id")]
    init <- data.frame(time = attr(object, "timeRange")[1L],
                       type = NA_character_, id = NA_character_,
                       nSusceptible = n, nInfectious = m, nRemoved = 0L)
    cumulatedReSusceptibility <- cumsum(counters[["type"]] == "S")
    cumulatedInfections <- cumsum(counters[["type"]] == "I")
    cumulatedRemovals <- cumsum(counters[["type"]] == "R")
    counters[["nSusceptible"]] <-
        init[["nSusceptible"]] - cumulatedInfections + cumulatedReSusceptibility
    counters[["nInfectious"]] <-
        init[["nInfectious"]]  + cumulatedInfections - cumulatedRemovals
    counters[["nRemoved"]] <-
        init[["nRemoved"]]     + cumulatedRemovals   - cumulatedReSusceptibility
    counters <- rbind(init, counters)
    rownames(counters) <- NULL

    ### return the components in a list
    res <- list(type = epitype, size = n - sum(isNeverInfected),
        initiallyInfected = initiallyInfected, neverInfected = neverInfected,
        coordinates = coordinates, byID = byID, counters = counters)
    class(res) <- "summary.epidata"
    attr(res, "eventTimes") <- attr(object, "eventTimes")
    attr(res, "timeRange") <- attr(object, "timeRange")
    res
}



################################################################################
# PRINT METHOD FOR 'EPIDATA' OBJECTS
################################################################################

print.epidata <- function (x, ...)
{
    cat("\nHistory of an epidemic\n")
    cat("Number of individuals:", nlevels(x[["id"]]), "\n")
    cat("Time range:", paste(attr(x, "timeRange"), collapse = " -- "), "\n")
    cat("Number of infections:", length(attr(x, "eventTimes")), "\n\n")
    print.data.frame(x, ...)
    cat("\n")
    invisible(x)
}


################################################################################
# PRINT METHOD FOR THE SUMMARY OF 'EPIDATA' OBJECTS
################################################################################

print.summary.epidata <- function(x, ...)
{
    cat("\nAN", x$type, "EPIDEMIC\n")
    cat("  Time range:", paste(attr(x, "timeRange"), collapse = " -- "), "\n")
    cat("  Number of individuals:", nlevels(x$initiallyInfected), "\n")
    cat(" ", length(x$initiallyInfected), "initially infected individuals")
    if (length(x$initiallyInfected) > 0) {
        cat(":\n    ")
        str(as.character(x$initiallyInfected), give.head = FALSE, vec.len = 100,
            strict.width = "wrap", indent.str = "  ")
    } else cat("\n")
    cat(" ", length(x$neverInfected), "never infected individuals")
    if (length(x$neverInfected) > 0) {
        cat(":\n    ")
        str(as.character(x$neverInfected), give.head = FALSE, vec.len = 100,
            strict.width = "wrap", indent.str = "  ")
    } else cat("\n")
    cat("  Size of the epidemic:", x$size, "\n")
    if (x$type %in% c("SIRS", "SIS")) {
        cat("  Number of infections:", length(attr(x, "eventTimes")), "\n")
    }
    dimc <- dim(x$counters)
    cat("\n$ counters ('data.frame',", dimc[1L], "x", dimc[2L], "):",
        "evolution of the epidemic:\n")
    counters2print <- if (dimc[1] > 6L) {
            tmp <- format.data.frame(x$counters[c(1:3,1,dimc[1]-(1:0)),],
                                     na.encode = FALSE)
            tmp[4,] <- c("[....]", "", "", "", "", "")
            rownames(tmp)[4] <- ""
            as.matrix(tmp)
        } else { x$counters }
    print(counters2print, quote = FALSE, right = TRUE, na.print = "")
    cat("\n")
    invisible(x)
}
