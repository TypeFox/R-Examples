################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### plot-method for "epidataCS" objects
###
### Copyright (C) 2009-2015 Sebastian Meyer
### $Revision: 1507 $
### $Date: 2015-11-04 01:09:33 +0100 (Mit, 04. Nov 2015) $
################################################################################


plot.epidataCS <- function (x, aggregate = c("time", "space"), subset,
                            by = type, ...)
{
    aggregate <- match.arg(aggregate)
    FUN <- paste("epidataCSplot", aggregate, sep = "_")
    do.call(FUN, args = list(x = quote(x), subset = substitute(subset),
                             by = substitute(by), ...))
}


### plot.epidataCS(x, aggregate = "time") -> number of cases over time
## in case t0.Date is specified, hist.Date() is used and breaks must set in ... (e.g. "months")

epidataCSplot_time <- function (x, subset, by = type,
    t0.Date = NULL, breaks = "stgrid", freq = TRUE,
    col = rainbow(nTypes), cumulative = list(), add = FALSE, mar = NULL,
    xlim = NULL, ylim = NULL, xlab = "Time", ylab = NULL, main = NULL,
    panel.first = abline(h=axTicks(2), lty=2, col="grey"),
    legend.types = list(), ...)
{
    timeRange <- with(x$stgrid, c(start[1L], stop[length(stop)]))
    ## subset event marks
    eventMarks <- if (missing(subset)) {
        marks.epidataCS(x, coords = FALSE)
    } else {
        do.call(base::subset, list(
            x = quote(marks.epidataCS(x, coords = FALSE)),
            subset = substitute(subset)
        ))
    }
    if (nrow(eventMarks) == 0L) stop("no events left after 'subset'")
    ## extract the data to plot
    by <- substitute(by)
    eventTimesTypes <- eventMarks[c("time", "type")]
    eventTimesTypes$type <- if (is.null(by)) { # disregard event types
        factor("all")
    } else { # stratification of counts (default is to stack bars by event type)
        as.factor(eval(by, envir = eventMarks))
    }
    typeNames <- levels(eventTimesTypes$type)
    nTypes <- length(typeNames)
    if (!freq && nTypes > 1L)
        warning("a stacked barplot of multiple event types only makes sense for 'freq=TRUE'")

    ## default breaks at stop times of stgrid
    if (identical(breaks, "stgrid")) {
        breaks <- c(timeRange[1L], unique.default(x$stgrid$stop))
        if (any(eventTimesTypes$time < timeRange[1L])) {
            message("Note: ignoring events of the pre-history (before \"stgrid\")")
            eventTimesTypes <- base::subset(eventTimesTypes, time >= timeRange[1L])
            if (nrow(eventTimesTypes) == 0L) stop("no events left to plot")
        }
    }

    ## calculate cumulative numbers if requested
    if (is.list(cumulative)) {
        csums <- tapply(eventTimesTypes$time, eventTimesTypes["type"],
                        function (t) cumsum(table(t)), simplify=FALSE)
        if (!is.null(cumulative[["offset"]])) {
            stopifnot(is.vector(cumulative$offset, mode="numeric"),
                      length(cumulative$offset) == nTypes)
            csums <- mapply(FUN="+", csums, cumulative$offset,
                            SIMPLIFY=FALSE, USE.NAMES=TRUE)
        }
        if (is.null(cumulative[["axis"]])) cumulative[["axis"]] <- TRUE
    }
    
    eventTimesTypes$type <- as.integer(eventTimesTypes$type)
    typesEffective <- sort(unique(eventTimesTypes$type))
    col <- rep_len(col, nTypes)
    
    if (!is.null(t0.Date)) {
        stopifnot(length(t0.Date) == 1L)
        t0.Date <- as.Date(t0.Date)
        t0 <- timeRange[1L]
        if (is.numeric(breaks) && length(breaks) > 1L)  # transform to Date
            breaks <- t0.Date + (breaks - t0)
        if (is.null(xlim))
            xlim <- t0.Date + (timeRange - t0)
        if (missing(xlab) && is.character(breaks))
            xlab <- paste0("Time (", breaks, ")")
        eventTimesTypes$time <- t0.Date + as.integer(eventTimesTypes$time - t0)
        ## we need integer dates here because otherwise, if the last event
        ## occurs on the last day of a month, year, etc. (depending on
        ## 'breaks') with a fractional date (e.g. as.Date("2009-12-31") + 0.5),
        ## then the automatic 'breaks' (e.g., breaks = "months") will not cover
        ## the data (in the example, it will only reach until
        ## as.Date("2009-12-31")). The following would fail:
        ## data("imdepi"); plot(imdepi, t0.Date = "2002-01-15", breaks = "months")
    }
    
    gethistdata <- function (breaks, types = seq_len(nTypes)) {
        times <- eventTimesTypes$time[eventTimesTypes$type %in% types]
        if (is.null(t0.Date)) {
            hist(times, breaks=breaks, plot=FALSE, warn.unused=FALSE, ...)
        } else {
            hist(times, breaks=breaks, plot=FALSE, ...)
            ## warn.unused=FALSE is hard-coded in hist.Date
        }
    }
    histdata <- gethistdata(breaks=breaks)
    if (!is.null(t0.Date)) {
        ## hist.Date() drops the Date class, but we need it for later re-use
        class(histdata$breaks) <- "Date"
    }
    
    ## establish the basic plot window
    if (!add) {
        if (is.null(xlim)) xlim <- timeRange
        if (is.null(ylim)) {
            ylim <- range(0, histdata[[if (freq) "counts" else "density"]])
        }
        if (is.null(ylab)) {
            ylab <- if (freq) "Number of cases" else "Density of cases"
        }
        
        if (is.null(mar)) {
            mar <- par("mar")
            if (is.list(cumulative) && cumulative$axis) mar[4L] <- mar[2L]
        }
        opar <- par(mar = mar); on.exit(par(opar))
        plot(x=xlim, y=ylim, xlab=xlab, ylab=ylab, main=main, type="n", bty="n")
        force(panel.first)
    }

    ## plot histogram (over all types)
    suppressWarnings( # about wrong AREAS if breaks are non-equidistant
        plot(histdata, freq = freq, add = TRUE, col = col[typesEffective[1L]], ...)
    )
    if (!add)  # doesn't work as expected when adding to plot with cumulative axis
        box()  # because white filling of bars might overdraw the inital box

    ## add type-specific sub-histograms
    for (typeIdx in seq_along(typesEffective)[-1L]) {
        .histdata <- gethistdata(
            breaks = histdata$breaks, # have to use same breaks
            types = typesEffective[typeIdx:length(typesEffective)]
        )
        suppressWarnings( # about wrong AREAS if breaks are non-equidistant
            plot(.histdata, freq = freq, add = TRUE,
                 col = col[typesEffective[typeIdx]], ...)
        )
    }

    ## optionally add cumulative number of cases
    if (is.list(cumulative)) {
        aT2 <- axTicks(2)
        div <- length(aT2) - 1L
        darken <- function (col, f = 0.6)
            apply(X = col2rgb(col, alpha = TRUE), MARGIN = 2L,
                  FUN = function (x) rgb(f*x[1L], f*x[2L], f*x[3L], x[4L],
                                         maxColorValue = 255))
        cumulative <- modifyList(
            list(maxat = ceiling(max(unlist(csums))/div)*div,
                 col = darken(col), lwd = 3, axis = TRUE,
                 lab = "Cumulative number of cases"),
            cumulative)
        csum2y <- function (x) x / cumulative$maxat * aT2[length(aT2)]
        for (typeIdx in typesEffective) {
            .times <- as.numeric(names(csums[[typeIdx]]))
            lines(if (is.null(t0.Date)) .times else t0.Date + .times - t0,
                  csum2y(csums[[typeIdx]]), lwd=cumulative$lwd,
                  col=cumulative$col[typeIdx])
        }
        if (cumulative$axis) {
            axis(4, at=aT2, labels=aT2/aT2[length(aT2)]*cumulative$maxat)
            mtext(cumulative$lab, side=4, line=3, las=0)
        }
    }

    ## optionally add legend
    if (is.list(legend.types) && length(typesEffective) > 1) {
        legend.types <- modifyList(
            list(x="topleft", legend=typeNames[typesEffective],
                 title=deparse(by, nlines = 1), fill=col[typesEffective]),
            legend.types)
        do.call("legend", legend.types)
    }
    
    invisible(histdata)
}


### plot.epidataCS(x, aggregate = "space") -> spatial point pattern

epidataCSplot_space <- function (x, subset, by = type, tiles = x$W, pop = NULL,
    cex.fun = sqrt, points.args = list(), add = FALSE,
    legend.types = list(), legend.counts = list(), sp.layout = NULL, ...)
{
    ## extract the points to plot
    events <- if (missing(subset)) {
        x$events
    } else { # calls sp:::subset.Spatial
        eval(substitute(base::subset(x$events, subset=.subset),
                        list(.subset=substitute(subset))))
    }
    ## should the plot distinguish between different event types?
    by <- substitute(by)
    events@data$type <- if (is.null(by)) { # disregard event types
        factor("all")
    } else { # default is to distinguish points by event type
        as.factor(eval(by, envir = events@data))
    }
    typeNames <- levels(events$type)
    nTypes <- length(typeNames)
    eventCoordsTypes <- data.frame(
        coordinates(events), type = as.integer(events$type),
        row.names = NULL, check.rows = FALSE, check.names = FALSE)

    ## count events by location and type
    eventCoordsTypesCounts <- if (is.null(pop)) {
        countunique(eventCoordsTypes)
    } else {
        ## work with "SpatialPolygons" -> spplot()
        events$COUNT <- multiplicity(eventCoordsTypes)
        events[!duplicated(eventCoordsTypes), c("type", "COUNT")]
    }
    pointCounts <- eventCoordsTypesCounts$COUNT
    countsLegend <- unique(round(10^(do.call("seq", c(
        as.list(log10(range(pointCounts))), list(length.out=5)
    )))))
    typesEffective <- sort(unique(eventCoordsTypesCounts$type))

    ## point style
    colTypes <- list(...)[["colTypes"]]  # backwards compatibility for < 1.8
    if (is.null(colTypes)) {
        colTypes <- rainbow(nTypes)
    } else warning("argument 'colTypes' is deprecated; ",
                   "use 'points.args$col' instead")
    points.args <- modifyList(list(pch=1, col=colTypes, lwd=1, cex=0.5),
                              points.args)
    styleArgs <- c("pch", "col", "lwd")
    points.args[styleArgs] <- lapply(points.args[styleArgs],
                                     rep_len, length.out=nTypes)

    ## select style parameters according to the events' types
    points.args_pointwise <- points.args
    points.args_pointwise[styleArgs] <- lapply(
        points.args_pointwise[styleArgs], "[",
        eventCoordsTypesCounts$type)
    points.args_pointwise$cex <- points.args_pointwise$cex * cex.fun(pointCounts)
    
    ## plot
    if (is.null(pop)) {
        ## classical plotting system
        if (!add) plot(tiles, ...)
        do.call("points", c(alist(x=eventCoordsTypesCounts[,1:2,drop=FALSE]),
                            points.args_pointwise))
        ## optionally add legends
        if (is.list(legend.types) && length(typesEffective) > 1) {
            legend.types <- modifyList(
                list(x="topright", legend=typeNames[typesEffective],
                     title=deparse(by, nlines = 1),
                     #pt.cex=points.args$cex, # better use par("cex")
                     pch=points.args$pch[typesEffective],
                     col=points.args$col[typesEffective],
                     pt.lwd=points.args$lwd[typesEffective]),
                legend.types)
            do.call("legend", legend.types)
        }
        if (is.list(legend.counts) && any(pointCounts > 1)) {
            if (!is.null(legend.counts[["counts"]])) {
                countsLegend <- as.vector(legend.counts[["counts"]], mode="integer")
                legend.counts[["counts"]] <- NULL
            }
            legend.counts <- modifyList(
                list(x="bottomright", bty="n", legend=countsLegend,
                     pt.cex=points.args$cex * cex.fun(countsLegend),
                     pch=points.args$pch[1L],
                     col=if(length(unique(points.args$col)) == 1L)
                             points.args$col[1L] else 1,
                     pt.lwd=points.args$lwd[1L]),
                legend.counts)
            do.call("legend", legend.counts)
        }
        invisible()
    } else {
        if (!is(tiles, "SpatialPolygonsDataFrame")) {
            stop("'pop' requires 'tiles' to be a \"SpatialPolygonsDataFrame\"")
        }
        ## grid plotting system -> spplot()
        layout.points <- c(list("sp.points", eventCoordsTypesCounts),
                           points.args_pointwise)
        ## optional legend definitions
        legend.types <- if (is.list(legend.types) && length(typesEffective) > 1) {
            legend.types <- modifyList(
                list(corner = c(1, 1), # "topright"
                     title = deparse(by, nlines = 1), cex.title = 1, border = TRUE,
                     points = list(
                         pch = points.args$pch[typesEffective],
                         col = points.args$col[typesEffective],
                         lwd = points.args$lwd[typesEffective]
                     ),
                     text = list(typeNames[typesEffective])),
                legend.types
            )
            corner.types <- legend.types$corner
            legend.types$corner <- NULL
            list(inside = list(fun = lattice::draw.key(legend.types), corner = corner.types))
        }
        legend.counts <- if (is.list(legend.counts) && any(pointCounts > 1)) {
            if (!is.null(legend.counts[["counts"]])) {
                countsLegend <- as.vector(legend.counts[["counts"]], mode="integer")
                legend.counts[["counts"]] <- NULL
            }
            legend.counts <- modifyList(
                list(corner = c(1,0), # "bottomright"
                     points = list(
                         cex = points.args$cex * cex.fun(countsLegend),
                         pch = points.args$pch[1L],
                         col = if(length(unique(points.args$col)) == 1L)
                             points.args$col[1L] else 1,
                         lwd = points.args$lwd[1L]
                     ),
                     text = list(as.character(countsLegend)),
                     padding.text=2, between=0),
                legend.counts
            )
            corner.counts <- legend.counts$corner
            legend.counts$corner <- NULL
            list(inside = list(fun = lattice::draw.key(legend.counts), corner = corner.counts))
        }
        ## create the plot
        spplot(obj = tiles, zcol = pop,
               sp.layout = c(list(layout.points), sp.layout),
               legend = c(legend.types, legend.counts), ...)
    }
}
