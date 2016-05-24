################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### The plot-method for "epidata" (via plot.summary.epidata) shows the evolution
### of the numbers of susceptible, infectious and recovered individuals.
### The extra function "stateplot" shows the event history of one individual.
###
### Copyright (C) 2008-2009, 2013-2014 Sebastian Meyer
### $Revision: 1080 $
### $Date: 2014-10-19 00:00:08 +0200 (Son, 19. Okt 2014) $
################################################################################

plot.epidata <- function(x, ...)
{
    sx <- summary(x)
    plot.summary.epidata(sx, ...)
}

plot.summary.epidata <- function (x,
    lty = c(2,1,3), lwd = 2,
    col = c("#1B9E77", "#D95F02", "#7570B3"), col.hor = col, col.vert = col,
    xlab = "Time", ylab = "Number of individuals", xlim = NULL, ylim = NULL,
    legend.opts = list(), do.axis4 = NULL, panel.first = grid(),
    rug.opts = list(), which.rug = c("infections", "removals",
    "susceptibility", "all"), ...)
{
    counters <- x[["counters"]]
    type <- x[["type"]]

    n <- counters[1L,"nSusceptible"]
    m <- counters[1L,"nInfectious"]
    N <- n + m
    times <- counters[-1L,"time"]
    if (missing(lty)) {
        lty <- c(2, 1, 3 * (type %in% c("SIR","SIRS")))
    }
    recycle3 <- function (xnam)
        assign(xnam, rep(get(xnam), length.out = 3), inherits = TRUE)
    for(varname in c("lty", "lwd", "col", "col.hor", "col.vert"))
        recycle3(varname)
    if (is.null(xlim)) {
        xlim <- attr(x, "timeRange")
        if (xlim[2] == Inf) xlim[2] <- times[length(times)]
    }
    if (is.null(ylim))
        ylim <- c(0, max(
            (lty[1] > 0) * {if (type %in% c("SIRS", "SIS")) N else n},
            (lty[2] > 0) * max(counters$nInfectious),
            (lty[3] > 0) * max(counters$nRemoved)
            ))

    # basic plotting frame
    plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab,
         panel.first = panel.first, ...)
    abline(h = c(0, N), col = "grey")

    # for real xlim in lines.stepfun (see 'dr' adjustment in plot.stepfun code)
    fakexlim <- c(1,2) * (xlim[2] + 2*xlim[1])/3 - c(0,xlim[1])
    # this isn't nice, a user argument 'dr' in plot.stepfun would be appreciated
    
    # add #Susceptibles
    if (all(counters$nSusceptible == n)) {
        lines(x = xlim, y = c(n,n),
              lty = lty[1], lwd = lwd[1], col = col.hor[1], ...)
    } else {
        lines(stepfun(times, counters$nSusceptible), xlim = fakexlim,
              lty = lty[1], lwd = lwd[1], col.hor = col.hor[1],
              col.vert = col.vert[1], do.points = FALSE, ...)
    }

    # add #Infected
    if (all(counters$nInfectious == m)) {
        lines(x = xlim, y = c(m,m),
              lty = lty[2], lwd = lwd[2], col = col.hor[2], ...)
    } else {
        lines(stepfun(times, counters$nInfectious),
              xlim = fakexlim, lty = lty[2], lwd = lwd[2], col.hor = col.hor[2],
              col.vert = col.vert[2], do.points = FALSE, ...)
    }

    # add #Removed
    if (all(counters$nRemoved == 0)) {
        lines(x = xlim, y = c(0,0),
              lty = lty[3], lwd = lwd[3], col = col.hor[3], ...)
    } else {
        lines(stepfun(times, counters$nRemoved),
              xlim = fakexlim, lty = lty[3], lwd = lwd[3], col.hor = col.hor[3],
              col.vert = col.vert[3], do.points = FALSE, ...)
    }

    # add special annotations
    if (is.null(do.axis4)) do.axis4 <- type == "SIR"
    if (do.axis4) {
        finalvalues <- counters[nrow(counters), c("nSusceptible", "nRemoved")]
        axis(4, at = finalvalues[lty[c(1,3)] > 0], font = 2, ...)
    }
    if (is.list(rug.opts)) {
        if (is.null(rug.opts$ticksize)) rug.opts$ticksize <- 0.02
        if (is.null(rug.opts$quiet)) rug.opts$quiet <- TRUE
        which.rug <- match.arg(which.rug)
        if (is.null(rug.opts$col)) rug.opts$col <-
            switch(which.rug, all = 1, infections = col.hor[2],
                   removals = col.hor[3], susceptibility = col.hor[1])
        rugLocations <- switch(which.rug,
            all = times, infections = attr(x, "eventTimes"),
            removals = counters$time[counters$type == "R"],
            susceptibility = counters$time[counters$type == "S"]
        )
        if (length(rugLocations) > 0) {
            do.call(rug, c(list(x = rugLocations), rug.opts))
        }
    }
    if (is.list(legend.opts)) {
        legend.opts <- modifyList(
            list(x = "topright", bty = "n", inset = c(0,0.02),
                 legend = c("susceptible", "infectious", "removed")[lty>0],
                 lty = lty[lty>0], lwd = lwd[lty>0], col = col.hor[lty>0]),
            legend.opts)
        do.call(legend, legend.opts)
    }
    invisible(as.matrix(
        counters[c("time", "nSusceptible", "nInfectious", "nRemoved")]
    ))
}


################################################################################
# PLOT THE STATE CHANGES OF ONE INDIVIDUAL OF "epidata"
# ... will be passed to the plot function (stepfun or curve),
# e.g. add, xlim, ylim, main, xlab, ylab, ...
################################################################################

stateplot <- function(x, id, ...)
{
    sx <- getSummary(x, class = "epidata")

    .id <- as.character(id)
    if (length(.id) != 1) {
        stop ("'id' must have length 1")
    }
    initiallyInfected <- sx[["initiallyInfected"]]
    if (! .id %in% levels(initiallyInfected)) {
        stop ("invalid 'id', does not exist in 'x'")
    }
    isInitiallyInfected <- .id %in% initiallyInfected
    
    counters <- sx[["counters"]]
    states <- levels(counters[["type"]])
    
    path <- counters[which(counters$id == .id), c("time", "type")]
    # remove pseudo-R-events, which come before S-event
    directSevents <- which(duplicated(path[["time"]]))
    path_noPseudoR <- if (length(directSevents)) {
            path[-(directSevents-1), ]
        } else {
            path
        }
    
    pathfunc <-
        if (nrow(path_noPseudoR) > 0) {
            stepfun(
                x = path_noPseudoR[["time"]],
                y = c(1+isInitiallyInfected, unclass(path_noPseudoR[["type"]])),
                right = FALSE
            )
        } else {
            function(t) rep(1+isInitiallyInfected, length(t))
        }
    
    # plot it
    dotargs <- list(...)
    nms <- names(dotargs)
    if(! "xlab" %in% nms) dotargs$xlab <- "time"
    if(! "ylab" %in% nms) dotargs$ylab <- "state"
    if(! "main" %in% nms) dotargs$main <- ""
    if(! "xlim" %in% nms) dotargs$xlim <- attr(sx, "timeRange")
    if(! "xaxs" %in% nms) dotargs$xaxs <- "i"
    if(! "do.points" %in% nms && inherits(pathfunc, "stepfun")) {
        dotargs$do.points <- FALSE
    }
    do.call("plot", args = c(list(x = pathfunc, yaxt = "n"), dotargs))
    axis(2, at = seq_along(states), labels = states)

    invisible(pathfunc)
}
