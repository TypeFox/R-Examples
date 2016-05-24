################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### animate-method for "epidataCS" objects
### It respects the ani.options() "interval" and "nmax" of the animation
### package, and it is advisable to use it within saveHTML() or similar
###
### Copyright (C) 2009-2014 Sebastian Meyer
### $Revision: 1096 $
### $Date: 2014-10-30 11:59:12 +0100 (Don, 30. Okt 2014) $
################################################################################


## three types:
## time.spacing=NULL: sequential snapshots at all event times
## time.spacing=scalar: snapshots with given time step (and timer)
## time.spacing=NA: time step is determined such that "nmax" snapshots result

animate.epidataCS <- function (object, interval = c(0,Inf), time.spacing = NULL,
    nmax = NULL, sleep = NULL, legend.opts = list(), timer.opts = list(),
    pch = 15:18, col.current = "red", col.I = "#C16E41", col.R = "#B3B3B3",
    col.influence = NULL, main = NULL, verbose = interactive(), ...)
{
    stopifnot(is.numeric(interval), length(interval) == 2L)
    with.animation <- requireNamespace("animation", quietly = TRUE)
    if (is.null(sleep)) {
        sleep <- if (with.animation) animation::ani.options("interval") else 0.1
        ## we cannot set this as default function argument, because we don't
        ## want to depend on package "animation" (surveillance only suggests it)
    }
    if (is.null(nmax)) {
        nmax <- if (with.animation) animation::ani.options("nmax") else Inf
    }
    s <- summary(object)
    removalTimes <- s$eventTimes + object$events$eps.t
    eventCoordsTypes <- cbind(s$eventCoords, type = s$eventTypes)
    pch <- rep_len(pch, s$nTypes)
    typeNames <- names(s$typeTable)
    multitype <- length(typeNames) > 1L

    # set default legend options
    doLegend <- if (is.list(legend.opts)) {
        if (is.null(legend.opts[["x"]])) legend.opts$x <- "topright"
        if (is.null(legend.opts$title))  legend.opts$title <-
            if (multitype) "type" else "state"
        if (is.null(legend.opts$legend)) { legend.opts$legend <-
            if (multitype) typeNames else c("infectious", if (!is.na(col.R)) "removed")
        }
        if (is.null(legend.opts$col)) { legend.opts$col <-
            if (multitype) col.current else c(col.I, if (!is.na(col.R)) col.R)
        }
        if (is.null(legend.opts$pch)) legend.opts$pch <- pch
        TRUE
    } else FALSE

    # set default timer options
    doTimer <- if (is.list(timer.opts)) {
        if (is.null(timer.opts[["x"]]))  timer.opts$x <- "bottomright"
        if (is.null(timer.opts$title))   timer.opts$title <- "time"
        if (is.null(timer.opts$box.lty)) timer.opts$box.lty <- 0
        if (is.null(timer.opts$adj))     timer.opts$adj <- c(0.5,0.5)
        if (is.null(timer.opts$inset))   timer.opts$inset <- 0.01
        if (is.null(timer.opts$bg))      timer.opts$bg <- "white"
        TRUE
    } else FALSE

    # wrapper for 'points' with specific 'cex' for multiplicity
    multpoints <- function (tableCoordsTypes, col) {
        tableMult <- countunique(tableCoordsTypes)
        points(tableMult[,1:2,drop=FALSE], pch = pch[tableMult[,"type"]],
               col = col, cex = sqrt(1.5*tableMult[,"COUNT"]/pi) * par("cex"))
    }
    # functions returning if events are in status I or R at time t
    I <- function (t) s$eventTimes <= t & removalTimes >= t
    R <- function (t) removalTimes < t

    sequential <- is.null(time.spacing)  # plot observed infections sequentially
    if (!sequential) stopifnot(length(time.spacing) == 1L)
    timeGrid <- if (sequential) unique(s$eventTimes) else {
        start <- max(s$timeRange[1], interval[1])
        end <- min(interval[2], s$timeRange[2],
            max(removalTimes) + if (is.na(time.spacing)) 0 else time.spacing)
        if (is.na(time.spacing)) {
            if (!is.finite(nmax)) {
                stop("with 'time.spacing=NA', 'nmax' must be finite")
            }
            seq(from = start, to = end, length.out = nmax)
        } else {
            tps <- seq(from = start, to = end, by = time.spacing)
            if (length(tps) > nmax) {
                message("Generating only the first ",
                        sQuote(if (with.animation) "ani.options(\"nmax\")" else "nmax"),
                        " (=", nmax, ") snapshots")
                head(tps, nmax)
            } else tps
        }
    }
    .info <- format.info(timeGrid)
    timerformat <- paste0("%", .info[1], ".", .info[2], "f")

    # animate
    loopIndex <- if (!sequential) timeGrid else {
        idxs <- which(s$eventTimes >= interval[1] & s$eventTimes <= interval[2])
        if (length(idxs) > nmax) {
            message("Generating only the first ",
                    sQuote(if (with.animation) "ani.options(\"nmax\")" else "nmax"),
                    " (=", nmax, ") events")
            head(idxs, nmax)
        } else idxs
    }
    told <- -Inf
    if (verbose)
        pb <- txtProgressBar(min=0, max=max(loopIndex), initial=0, style=3)
    for(it in loopIndex) {
        t <- if (sequential) s$eventTimes[it] else it
        infectious <- I(t)
        removed <- R(t)
        plot(object$W, ...)             # FIXME: use default lwd = 2
        title(main = main)
        if (doLegend) do.call(legend, legend.opts)
        if (doTimer) {
            ttxt <- sprintf(timerformat, t)
            do.call(legend, c(list(legend = ttxt), timer.opts))
        }
        if (!is.null(col.influence)) {
            iRids <- which(infectious)
            if (sequential) setdiff(iRids, it)
            for(j in iRids) {
                iR <- shift.owin(object$events@data$.influenceRegion[[j]],
                                 s$eventCoords[j,])
                plot(iR, add = TRUE, col = col.influence, border = NA)
            }
        }
        rTable <- eventCoordsTypes[removed,,drop=FALSE]
        if (nrow(rTable) > 0L) multpoints(rTable, col = col.R)
        iTable <- eventCoordsTypes[infectious,,drop=FALSE]
        if (nrow(iTable) > 0L) multpoints(iTable, col = col.I)
        infectiousNew <- if (sequential) it else infectious & !I(told)
        iTableNew <- eventCoordsTypes[infectiousNew,,drop=FALSE]
        if (nrow(iTableNew) > 0L) multpoints(iTableNew, col = col.current)
        told <- t
        if (verbose) setTxtProgressBar(pb, it)
        if (dev.interactive()) Sys.sleep(sleep)
    }
    if (verbose) close(pb)
    ## if (dev.interactive())
    ##     message("Note: use facilities of the \"animation\" package, e.g.,\n",
    ##             "      saveHTML() to view the animation in a web browser.")
    invisible(NULL)
}
