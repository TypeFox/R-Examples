################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Spatial and temporal tie-breaking of events
###
### Copyright (C) 2012-2014 Sebastian Meyer
### $Revision: 1005 $
### $Date: 2014-09-10 23:55:11 +0200 (Mit, 10. Sep 2014) $
################################################################################


## epidataCS-method
## makes use of untie.default (in time) and untie.matrix (in space)
untie.epidataCS <- function (x,
                             amount = list(t=NULL, s=NULL),
                             minsep = list(t=0, s=0),
                             direction = "left", keep.sources = FALSE,
                             ..., verbose = FALSE)
{
    stopifnot(is.list(amount), !is.null(names(amount)),
              is.list(minsep), !is.null(names(minsep)))
    minsep <- modifyList(list(t=0, s=0), minsep)
    do.spatial <- pmatch("s", names(amount), nomatch=0L) > 0L
    do.temporal <- pmatch("t", names(amount), nomatch=0L) > 0L
    if (!do.spatial && !do.temporal) {
        stop("no amounts specified, nothing to do")
    }

    ## Generate new events data frame
    events <- marks.epidataCS(x, coords=FALSE)
    newcoords <- if (do.spatial) {      # untie spatial coordinates
        untie.matrix(coordinates(x$events), amount$s, minsep$s, constraint=x$W, ...)
    } else coordinates(x$events)
    if (do.temporal) {                  # untie event times
        ## by default, we shift event times (non-symmetrically) to the left such
        ## that the shifted versions potentially stay in the same BLOCK of
        ## endemic covariates (the CIF is left-continuous).
        events$time <- untie.default(events$time, amount$t, minsep$t,
                                     direction=direction, sort=TRUE, ...)
        ## FIXME: Does sort=TRUE always make sense?
        ##        maybe only sort in untie.default if amount < minsep?
    }

    ## Generate epidataCS object with new events
    coordinates(events) <- newcoords    # -> SpatialPointsDataFrame
    #proj4string(events) <- proj4string(x$W)  # "proj4string<-" might change the
                                        # string e.g. add +towgs84=0,0,0,0,0,0,0
    events@proj4string <- x$W@proj4string
    npoly <- attr(x$events$.influenceRegion, "nCircle2Poly")
    clipper <- attr(x$events$.influenceRegion, "clipper")
    if (is.null(clipper))  # epidataCS < 1.8-1
        clipper <- "polyclip"
    res <- as.epidataCS(events=events, stgrid=x$stgrid[,-1L], W=x$W,
                        qmatrix=x$qmatrix, nCircle2Poly=npoly,
                        clipper=clipper, verbose=verbose)
    if (keep.sources) {
        res$events$.sources <- x$events$.sources
    }

    ## Done
    res
}

## untie event times by uniform jittering
untie.default <- function (x, amount = NULL, minsep = 0,
                           direction = c("symmetric", "left", "right"),
                           sort = NULL, giveup = 1000, ...)
{
    stopifnot(is.numeric(x), is.vector(x))
    distx <- dist(x)
    isPosDist <- distx > 0
    if (all(isPosDist)) return(x)       # no ties
    direction <- match.arg(direction)
    if (is.null(sort))                  # sort if x was sorted
        sort <- identical(order(x, decreasing=FALSE), seq_along(x))

    if (any(isPosDist)) {
        minsepx <- min(distx[isPosDist])     # smallest positive distance
        amount.bound <- if (direction=="symmetric") minsepx/2 else minsepx
        if (is.null(amount)) {
            amount <- amount.bound
        } else if (sort && abs(amount) > amount.bound) {
            warning("'amount' should not be greater than ",
                    if (direction=="symmetric") "half of ",
                    "the minimum separation (", format(amount.bound), ")")
        }
    } else if (is.null(amount)) {
        stop("default 'amount' does not work with completely tied 'x'")
    }

    shiftFUN <- switch(direction,
        symmetric = function (x) x + runif(length(x), -amount, amount),
        right = function (x) x + runif(length(x), 0, amount),
        left = function (x) x - runif(length(x), 0, amount))
    res <- .untie(x, shiftFUN, minsep)
    
    if (sort) base::sort(res) else res
}

## untie spatial coordinates by moving them by vectors drawn uniformly from a
## disc of radius 'amount', optionally respecting a region (constraint)
## inside which the jittered points should be located (of course, the initial
## points must also obey this constraint), and a minimum separation 'minsep'
untie.matrix <- function (x, amount = NULL, minsep = 0,
                          constraint = NULL, giveup = 1000, ...)
{
    stopifnot(is.numeric(x), is.matrix(x))
    dimx <- dim(x)
    if (dimx[2L] <= 1L) {
        untie.default(c(x), amount, minsep, giveup=giveup)
    } else if (dimx[2L] > 2L) {
        stop("spatial tie-breaking is only implemented for 2D coordinates")
    }
    if (is.null(amount)) {
        distx <- dist(x)
        isPosDist <- distx > 0
        ## take half of smallest distance, which guarantees that new points
        ## will be closer to previously tied points than to others
        if (any(isPosDist)) amount <- min(distx[isPosDist]) / 2 else
        stop("default 'amount' does not work with a single location only")
    }
    if (!is.null(constraint)) {
        stopifnot(inherits(constraint, "SpatialPolygons"))
        proj4string(constraint) <- CRS(NA_character_)
        outOfConstraint <- function (x) {
            is.na(over(SpatialPoints(x), constraint))
        }
        if (any(outOfConstraint(x)))
            stop("some points of the matrix 'x' don't respect the 'constraint'")
    } else outOfConstraint <- NULL

    shiftFUN <- function (x) x + runifdisc(nrow(x), amount)
    .untie(x, shiftFUN, minsep, outOfConstraint, giveup=giveup)
}

## workhorse for both vector and matrix 'x'
.untie <- function (x, shiftFUN, minsep = 0, outOfConstraintFUN = NULL, giveup = 1000)
{
    x <- res <- as.matrix(x)
    move <- rep.int(TRUE, nrow(x))      # initially move _all_ points
    ntry <- 0L
    updateMoveExpr <- .updateMoveExpr(minsep>0, is.function(outOfConstraintFUN))
    while((nleft <- sum(move)) > 0L && ntry < giveup) {
        res[move,] <- shiftFUN(x[move,,drop=FALSE])
        ## determine for the moved points if they are too close to another point
        ## or fall out of constraint -> try again
        eval(updateMoveExpr)
        ntry <- ntry + 1L
    }
    if (ntry >= giveup)
        warning("could not obey 'constraint' and/or 'minsep' for some points")
    if (ncol(res) == 1) res[,1] else res
}

## check if points with index 'idx' are too close (< minsep) to any other points
## (this function could probably be made more efficient, especially for
## length(idx) << nrow(pts), since we actually don't need all pairwise distances
## calculated by dist() but only those related to the idx-points)
.tooClose <- function (pts, idx, minsep) {
    distpts <- as.matrix(dist(pts))
    diag(distpts) <- Inf
    rowSums(distpts[idx,,drop=FALSE] < minsep) > 0
}

## generate expression which updates logical vector 'move' (points left to move)
.updateMoveExpr <- function(doClose = FALSE, doConstraint = FALSE)
{
    if (!doClose && !doConstraint) return(expression(move[move] <- FALSE))
    exprClose <- expression(movedTooClose <- .tooClose(res, move, minsep))
    exprConstraint <- if (doClose) { # only need to check those not too close
        expression(
            movedOutOfConstraint <- rep.int(FALSE, nleft),
            if (any(!movedTooClose)) movedOutOfConstraint[!movedTooClose] <-
            outOfConstraintFUN(res[move,,drop=FALSE][!movedTooClose,,drop=FALSE])
            )
    } else {
        expression(
            movedOutOfConstraint <- outOfConstraintFUN(res[move,,drop=FALSE])
            )
    }
    c(if (doClose) exprClose, if (doConstraint) exprConstraint,
      switch(doClose + 2*doConstraint,
             expression(move[move] <- movedTooClose),
             expression(move[move] <- movedOutOfConstraint),
             expression(move[move] <- movedTooClose | movedOutOfConstraint)
             )
      )
}
