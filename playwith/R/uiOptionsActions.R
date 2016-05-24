## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer


initOptionsActions <- function(playState)
{
    hasArgs <- playState$accepts.arguments
    ## sensible default for time.mode based on data
    ## DISABLED -- do not want to evaluate data by default
    ## and anyway, it is easier now to zoom along x axis
    #if (hasArgs && playState$.args$missing_time.mode) {
    if (is.na(playState$time.mode)) {
        if (hasArgs) {
            ## detect default for time.mode based on data
            dat <- try(xyData(playState, space="packet 1"), silent = TRUE)
            playState$time.mode <-
                (is.list(dat) &&
                 (inherits(dat$x, "ts") ||
                  inherits(dat$x, "zoo") ||
                  is.somesortoftime(dat$x))
                 )
            ## once only ## TODO: better to wait for manual switch
            playState$.args$missing_time.mode <- FALSE
        } else {
            playState$time.mode <- FALSE
        }
    }
    with (playState$widgets, {
        timeScrollbar["sensitive"] <- playState$time.mode
        timeEntry["sensitive"] <- playState$time.mode
    })
    if (playState$time.mode)
        time.mode_init(playState)
}

updateOptionsActions <- function(playState)
{
    aGroup <- playState$actionGroups[["PlotActions"]]
    hasArgs <- playState$accepts.arguments
    ## Time Mode
    hasTimeVec <- !is.null(playState$time.vector)
    aGroup$getAction("TimeMode")$setSensitive(hasArgs || hasTimeVec)
    needInit <- (aGroup$getAction("TimeMode")$getActive() !=
                 playState$time.mode)
    aGroup$getAction("TimeMode")$setActive(isTRUE(playState$time.mode))
    if (needInit)
        time.mode_init(playState)
    ## update scrollbar etc
    time.mode_update(playState)
    ## global options
    aGroup <- playState$actionGroups[["GlobalActions"]]
    ## Annotations options
    aGroup$getAction("ClipAnnot")$setActive(isTRUE(playState$clip.annotations))
    aGroup$getAction("PageAnnot")$setActive(isTRUE(playState$page.annotation))
    ## Statusbar, toolbars, tooltips
    aGroup$getAction("ShowStatusbar")$setActive(isTRUE(playState$show.statusbar))
    aGroup$getAction("ShowToolbars")$setActive(isTRUE(playState$show.toolbars))
    aGroup$getAction("ShowTooltips")$setActive(isTRUE(playState$show.tooltips))
}

time.mode_handler <- function(widget, playState)
{
    playState$time.mode <- widget["active"]
    time.mode_init(playState)
    ## update scrollbar etc
    time.mode_update(playState)
}

time.mode_init <- function(playState)
{
    blockRedraws(with (playState$widgets, {
        timeScrollBox["visible"] <- TRUE
        timeScrollbar["sensitive"] <- playState$time.mode
        timeEntry["sensitive"] <- playState$time.mode
    }))
    ## store data range and class
    if (playState$time.mode) {
        if (is.null(playState$time.vector)) {
            xy <- xyData(playState, space="page")
            playState$tmp$time.mode.x.range <- range(as.numeric(xy$x))
            playState$tmp$time.mode.x.attr <- attributes(xy$x)
        }
    }
}

time.mode_update <- function(playState)
{
    if (!isTRUE(playState$time.mode)) return()
    #if (widget["active"] == FALSE) {
    #    widget["active"] <- TRUE ## triggers update
    #    return()
    #}
    blockRedraws({
        widg <- playState$widgets
        if (!is.null(playState$time.vector)) {
            x.pos <- playState$env$cur.index
            x.max <- length(playState$time.vector)
            x.jump <- playState$time.mode.page.incr
            if (is.null(x.jump)) x.jump <- round(log2(x.max))
            cur.time <- playState$env$cur.time
            if (inherits(cur.time, "yearqtr"))
              cur.time <- as.yearmon(cur.time)
            widg$timeEntry["text"] <- toString(cur.time)
            widg$timeScrollbar["adjustment"] <-
                gtkAdjustment(value=x.pos, lower=1, upper=x.max+1,
                              step.incr=1, page.incr=x.jump, page.size=1)
            widg$timeScrollbar$setValue(x.pos) ## need this (bug?)
            return()
        }
        x.range <- playState$tmp$time.mode.x.range
        x.lim <- rawXLim(playState)
        x.page <- abs(diff(x.lim))
        x.page <- min(x.page, abs(diff(x.range)))
        x.pos <- min(x.lim)
        x.pos <- max(x.pos, min(x.range))
        x.pos <- min(x.pos, max(x.range))
        ## format x limits for text box
        xlim <- signif(x.lim, 6)
        class(x.lim) <- playState$tmp$time.mode.x.attr$class
        if ("yearmon" %in% class(x.lim))
            x.lim <- as.yearmon(x.lim) ## to round correctly
        if ("yearqtr" %in% class(x.lim))
            x.lim <- as.yearmon(x.lim) ## to format as "%b %Y"
        if ("POSIXt" %in% class(x.lim))
            attr(x.lim, "tzone") <- playState$tmp$time.mode.x.attr$tzone
        if ("factor" %in% class(x.lim))
            attr(x.lim, "levels") <- playState$tmp$time.mode.x.attr$levels
        widg$timeEntry["text"] <- paste(format(x.lim), collapse=" to ")
        ## set up scrollbar
        widg$timeScrollbar["adjustment"] <-
            gtkAdjustment(value=x.pos, lower=min(x.range), upper=max(x.range),
                          step.incr=x.page/2, page.incr=x.page, page.size=x.page)
        widg$timeScrollbar$setValue(x.pos) ## need this (bug?)
    })
}

time.mode_scrollbar_handler <- function(widget, playState)
{
    if (!isTRUE(playState$tmp$plot.ready)) return()
    newLim <- widget$getValue()
    if (!is.null(playState$time.vector)) {
        newLim <- round(newLim)
        playState$env$cur.index <- newLim
        playState$env$cur.time <- playState$time.vector[newLim]
        playReplot(playState)
        return()
    }
    newLim[2] <- newLim + widget["adjustment"]["page-size"]
    if (widget["adjustment"]["page-size"] == 0) return() ## sanitycheck
                                        #oldLim <- rawXLim(playState)
                                        #if (min(oldLim) == min(newLim)) return()
    newLim <- round(newLim, 7)
    if (any(!is.finite(newLim))) return()
    rawXLim(playState) <- newLim
    playReplot(playState)
}

time.mode_entry_handler <- function(widget, playState)
{
    if (!isTRUE(playState$tmp$plot.ready)) return()
    if (!is.null(playState$time.vector)) {
        newLim <- widget["text"]
        time.vector <- playState$time.vector
        max.x <- length(time.vector)
        cls <- class(time.vector)
        if ("POSIXt" %in% cls) newLim <- try(as.POSIXct(newLim))
        else if ("Date" %in% cls) newLim <- try(as.Date(newLim))
        else if ("yearmon" %in% cls) newLim <- try(as.yearmon(newLim, "%b %Y"))
        else if ("yearqtr" %in% cls) newLim <- try(as.yearqtr(as.yearmon(newLim, "%b %Y")))
        if (inherits(newLim, "try-error")) {
            ## treat it as an index into time.vector
            cur.index <- try(as.integer(widget["text"]), silent=TRUE)
            if (inherits(cur.index, "try-error")) {
                ## give up
                gmessage.error(conditionMessage(newLim))
                return()
            }
        } else {
            newLim <- as.numeric(newLim)
            cur.index <- findInterval(newLim, time.vector)
        }
        cur.index <- max(1, min(max.x, cur.index))
        playState$env$cur.index <- cur.index
        playState$env$cur.time <- time.vector[cur.index]
        playReplot(playState)
        return()
    }
    newLim <- strsplit(widget["text"], " to ")[[1]]
    if ((length(newLim) != 2)) {
        gmessage.error("Give bounds in form \"LOWER to UPPER\".")
        return()
    }
    x.attr <- playState$tmp$time.mode.x.attr
    cls <- x.attr$class
    if ("POSIXt" %in% cls) newLim <- as.POSIXct(newLim)
    else if ("Date" %in% cls) newLim <- as.Date(newLim)
    else if ("yearmon" %in% cls) newLim <- as.yearmon(newLim, "%b %Y")
    else if ("yearqtr" %in% cls) newLim <- as.yearqtr(as.yearmon(newLim, "%b %Y"))
    else if ("integer" %in% cls) newLim <- as.integer(newLim)
    newLim <- as.numeric(newLim)
    if (any(!is.finite(newLim))) return()
    rawXLim(playState) <- newLim
    playReplot(playState)
}

clip.annotations_handler <- function(widget, playState)
{
    playState$clip.annotations <- widget["active"]
}

page.annotation_handler <- function(widget, playState)
{
    playState$page.annotation <- widget["active"]
}

show.statusbar_handler <- function(widget, playState) {
    playState$widgets$statusbarBox["visible"] <- widget["active"]
}

show.toolbars_handler <- function(widget, playState) {
    playState$widgets$leftToolbar["visible"] <- widget["active"]
}

show.tooltips_handler <- function(widget, playState)
{
    playState$show.tooltips <- widget["active"]
    ## start it (or stop it)
    gtkmain_handler(playState=playState)
}
