################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Two types of spatio-temporal animations of "epidata" are supported:
### - sequential plots regardless of time between events (i.e. only ordering)
### - chronological animation with timer
###
### Copyright (C) 2008-2009, 2012, 2014 Sebastian Meyer
### $Revision: 1096 $
### $Date: 2014-10-30 11:59:12 +0100 (Don, 30. Okt 2014) $
################################################################################

animate.epidata <- function (object, ...)
{
    s <- summary(object)
    animate.summary.epidata(s, ...)
}

animate.summary.epidata <- function (object,
    main = "An animation of the epidemic",
    pch = 19, col = c(3, 2, gray(0.6)), time.spacing = NULL,
    sleep = quote(5/.nTimes), legend.opts = list(), timer.opts = list(),
    end = NULL, generate.snapshots = NULL, ...)
{
    counters <- object[["counters"]]
    # remove pseudo-R-events, which come before S-event
    directSevents <- which(duplicated(counters[["time"]]))
    counters_noPseudoR <- if (length(directSevents)) {
            counters[-(directSevents-1), ]
        } else {
            counters
        }
    # remove initial row and keep essential columns
    eventTable <- counters_noPseudoR[-1, c("time", "type", "id")]
    eventTable[["type"]] <- unclass(eventTable[["type"]])  # get integer codes
    .nTimes <- nrow(eventTable)
    
    # extract initial individual information (id, at-risk, coordinates)
    coords <- object[["coordinates"]]
    d <- ncol(coords)
    if (d > 2L) {
        stop("spatial plotting in more than two dimensions is not implemented")
    } else if (d == 1L) {
        coords <- cbind(coords, 0)
    } else if (d == 0L) {
        stop ("'object' does not contain any defined coordinates")
    }
    
    # plot the initial state
    pch <- rep(pch, length.out = 3)
    col <- rep(col, length.out = 3)
    isInitiallyInfected <- rownames(coords) %in% object[["initiallyInfected"]]
    plot(coords, pch = ifelse(isInitiallyInfected, pch[2L], pch[1L]), 
                 col = ifelse(isInitiallyInfected, col[2L], col[1L]),
                 main = main, ...)
    if (is.list(legend.opts)) {
        if (is.null(legend.opts[["x",exact=TRUE]]))
            legend.opts$x <- "topright"
        if (is.null(legend.opts$legend))
            legend.opts$legend <- c("susceptible", "infectious", "removed")
        if (is.null(legend.opts$col)) legend.opts$col <- col
        if (is.null(legend.opts$pch)) legend.opts$pch <- pch
        do.call(legend, legend.opts)
    }
    
    # animate the epidemic by iteratively re-drawing points at the coordinates
    sleep <- eval(sleep)
    if (is.null(time.spacing)) { # plot events sequentially
        for(i in seq_len(.nTimes)) {
            if (dev.interactive()) Sys.sleep(sleep)
            tmp <- eventTable[i,]  # c(time, type, id)
            points(coords[as.character(tmp[["id"]]),,drop=FALSE],
                   pch = pch[tmp[["type"]]], col = col[tmp[["type"]]])
        }
    } else { # plot events chronologically
        if (is.null(end))
            end <- eventTable[.nTimes, "time"] + time.spacing
        timeGrid <- seq(from = time.spacing, to = end, by = time.spacing)
        timeWidth <- nchar(timeGrid[length(timeGrid)])
        timeDigits <- nchar(strsplit(as.character(time.spacing), ".",
            fixed = TRUE)[[1L]][2L])
        form <- paste("%", timeWidth, ".", timeDigits, "f", sep = "")
        if (is.list(timer.opts)) {
            if (is.null(timer.opts[["x",exact=TRUE]]))
                timer.opts$x <- "bottomright"
            if (is.null(timer.opts$title))   timer.opts$title <- "time"
            if (is.null(timer.opts$box.lty)) timer.opts$box.lty <- 0
            if (is.null(timer.opts$adj))     timer.opts$adj <- c(0.5,0.5)
            if (is.null(timer.opts$inset))   timer.opts$inset <- 0.01
            if (is.null(timer.opts$bg))      timer.opts$bg <- "white"
            do.call(legend, c(list(legend = sprintf(form, 0)), timer.opts))
        }
        oldtp <- tp <- attr(object, "timeRange")[1L]
        i <- 1L                   # to be used in the file argument in dev.print
        if (is.vector(generate.snapshots, mode="character") &&
            length(generate.snapshots) == 1L && requireNamespace("animation")) {
            img.name <- generate.snapshots
            ani.dev <- animation::ani.options("ani.dev")
            if (is.character(ani.dev)) ani.dev <- get(ani.dev)
            imgdir <- animation::ani.options("imgdir")
            imgtype <- animation::ani.options("ani.type")
            generate.snapshots <- list(
                device = ani.dev,
                file = quote(file.path(imgdir, paste0(img.name,i,".",imgtype))),
                width = animation::ani.options("ani.width"),
                height = animation::ani.options("ani.height")
            )
        }
        if (is.list(generate.snapshots)) {
            do.call(dev.print, generate.snapshots)
        }
        for(i in 1L+seq_along(timeGrid)) {
            tp <- timeGrid[i-1L]
            if (dev.interactive()) Sys.sleep(sleep)
            timeIndex <- which(eventTable[["time"]] > oldtp & eventTable[["time"]] <= tp)
            if (length(timeIndex) > 0L) {
                tmp <- eventTable[timeIndex,]  # c(time, type, id)
                points(coords[as.character(tmp[["id"]]),,drop=FALSE],
                       pch = pch[tmp[["type"]]], col = col[tmp[["type"]]])
            }
            if (is.list(timer.opts)) {
                do.call(legend, c(list(legend = sprintf(form,tp)), timer.opts))
            }
            oldtp <- tp
            if (is.list(generate.snapshots)) {
                do.call(dev.print, generate.snapshots)
            }
        }
    }
    invisible(NULL)
}
