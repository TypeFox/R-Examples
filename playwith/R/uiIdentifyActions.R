## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer


initIdentifyActions <- function(playState)
{
    guessLabels(playState)
}

updateIdentifyActions <- function(playState)
{
    ## draw persistent labels and brushed points
    canIdent <- playState$tmp$identify.ok
    if (canIdent) {
        drawLabels(playState)
        drawLinkedLocal(playState)
    }
    updateIdentifyActionStates(playState)
}

updateIdentifyActionStates <- function(playState)
{
    aGroup <- playState$actionGroups[["PlotActions"]]
    ## Identify etc
    canIdent <- playState$tmp$identify.ok
    aGroup$getAction("Identify")$setSensitive(canIdent)
    aGroup$getAction("IdTable")$setSensitive(canIdent)
    aGroup$getAction("FindLabels")$setSensitive(FALSE) #canIdent)
    hasIDs <- (length(playGetIDs(playState)) > 0)
    aGroup$getAction("SaveIDs")$setSensitive(hasIDs)
    ## Brush
    aGroup$getAction("Brush")$setSensitive(canIdent)
}

identifyCore <- function(playState, foo)
{
    if (!isTRUE(playState$tmp$identify.ok)) return()
    if (is.null(playState$labels)) return()
    if (is.null(foo$coords)) return()
    space <- foo$space
    data <- xyCoords(playState, space=foo$space)
    ## convert to log scale if necessary
    data <- dataCoordsToSpaceCoords(playState, data)
    if (length(data$x) == 0) return(FALSE)
    if (length(data$y) == 0) return(FALSE)
    coords <- foo$coords
    if (foo$is.click) {
        x <- coords$x[1]
        y <- coords$y[1]
        ppxy <- playDo(playState,
                       list(lx=convertX(unit(x, "native"), "points", TRUE),
                            ly=convertY(unit(y, "native"), "points", TRUE),
                            px=convertX(unit(data$x, "native"), "points", TRUE),
                            py=convertY(unit(data$y, "native"), "points", TRUE)),
                       space=foo$space)
        pdists <- with(ppxy, sqrt((px - lx)^2 + (py - ly)^2))
        ## all data points within 11 points
        which <- which(pdists < 11)
        if (length(which) == 0) return()
        ## order by distance from click
        which <- which[order(pdists[which])]
        ## account for multiple points (matrix values of data$x, data$y)
        n <- NROW(data$x)
        which <- unique(which %% n)
        which[which == 0] <- n

        idMenu <- gtkMenu()
        idMenu$popup(button=0, activate.time=gtkGetCurrentEventTime())
        item <- gtkMenuItem("Add label to plot:")
        item["sensitive"] <- FALSE
        idMenu$append(item)
        for (w in which) {
            datx <- data$x[[w]]
            daty <- data$y[[w]]
            pos <- with(ppxy, getTextPosition(x = lx - px[w],
                                              y = ly - py[w]))
            ss <- data$subscripts[[w]]
            if (is.null(ss)) ss <- w
            label <- toString(playState$labels[[ss]])
            item <- gtkMenuItem(label)
            idMenu$append(item)
            gSignalConnect(item, "activate",
                           function(widget, user.data) {
                               ss <- user.data$ss
                               pos <- user.data$pos
                               ## store newly identified points in playState
                               playSetIDs(playState, ss, pos = pos,
                                          type = "labelled",
                                          space = space,
                                          add = TRUE)
                           }, data = list(ss = ss, pos = pos))
        }
        idMenu$append(gtkSeparatorMenuItem())
        item <- gtkMenuItem("Right-click on point for detail")
        item["sensitive"] <- FALSE
        idMenu$append(item)
        ## show the menu
        while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
    } else {
        ## drag
        foo <- playSelectData(playState, foo = foo)
        if (is.null(foo)) return()
        if (length(foo$which) == 0) return()
        with(foo, {
            if (!is.click) pos <- 1
            ## store newly identified points in playState
            playSetIDs(playState, subscripts, pos = pos,
                       type = "labelled",
                       space = space,
                       add = TRUE)
        })
    }
}

id.table_handler <- function(widget, playState)
{
    dat <- getDataArg(playState)
    if (is.null(dat))
        dat <- as.data.frame(xyData(playState, space = "page")[c("x","y")])
    ## add row numbers and names
    dat <- cbind(row = 1:NROW(dat), names = rownames(dat),
                 dat)
    w <- gbasicdialog("Select cases to be brushed",
                 #widget = tabW,
                 handler = function(h, ...) {
                     ids <- svalue(tabW)#, index = TRUE)
                     dispose(h$obj)
                     playSetIDs(playState, ids)
                 })
    size(w) <- c(600, 400)
    tabW <- gtable(dat, multiple = TRUE, container = w)
    cur.ids <- playGetIDs(playState)
    if (length(cur.ids) > 0)
        svalue(tabW, index = TRUE) <- cur.ids
    visible(w, set = TRUE)
}

id.find_handler <- function(widget, playState) {
    ## TODO
    gmessage.error("not yet implemented")
}

set.labels_handler <- function(widget, playState)
{
    playFreezeGUI(playState)
    on.exit(playThawGUI(playState))
    ## widgets
    box <- ggroup(horizontal = FALSE)
    datArg <- getDataArg(playState, eval = FALSE)
    dat <- try(eval(datArg, playState$env))
    labcode <- NULL
    labdesc <- NULL
    if (!is.null(dat)) {
        labcode <- colnames(dat)
        labdesc <- colnames(dat)
    }
    rnCode <- deparseOneLine(call("rownames", datArg))
    labcode <- c(labcode,
                 "xyData()$x",
                 "xyData()$y",
                 'with(xyData(), paste(y, x, sep="@"))',
                 "NULL",
                 rnCode)
    labdesc <- c(labdesc,
                 "data x values",
                 "data y values",
                 "data y@x values",
                 "data subscripts",
                 rnCode)
    labradio <- gradio(labdesc, selected = length(labcode), container = box,
                       handler = function(h, ...) {
                           idx <- max(1, svalue(labradio, index=TRUE))
                           svalue(labedit) <- labcode[idx]
                       })
    labedit <- gedit(rnCode, container = box)
    ## show dialog
    gbasicdialog("Set labels to...", widget = box,
                 handler = function(h, ...) {
                     playDevSet(playState)
                     expr <- parse(text=svalue(labedit))
                     tmp <- tryCatch(
                             eval(expr, dat, playState$env),
                                     error=function(e)e)
                     ## check whether there was an error
                     if (inherits(tmp, "error")) {
                         gmessage.error(conditionMessage(tmp))
                     } else {
                         ## set labels
                         playState$labels <- tmp
                         ## and set default for playNewPlot
                         playState$.args$labels <- tmp
                         ## redraw if needed
                         if (length(playGetIDs(playState, type = "labelled")))
                             playReplot(playState)
                     }
                 })
}

save.ids_handler <- function(widget, playState)
{
    name <- ginput(paste("Save subscripts of labelled /",
                         "brushed points to variable:"),
                   title = "Save IDs", text = "myIds")
    if ((length(name) == 0) || (nchar(name) == 0))
        return()
    playDevSet(playState)
    assign(name, playGetIDs(playState), globalenv())
}

brushCore <- function(playState, foo, add = TRUE)
{
    foo <- playSelectData(playState, foo = foo)
    if (is.null(foo)) return()
    if (length(foo$which) == 0) {
        ## no data points selected
        if (add == FALSE)
            playClear(playState, type = "brushed")
    } else {
        playSetIDs(playState, foo$subscripts,
                   add = add)
    }
}

