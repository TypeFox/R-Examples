# playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

playSelectData <-
    function(playState = playDevCur(),
             prompt = paste(
             "Click or drag to select data points;",
             "Right-click or Esc to cancel."),
             scales = "dynamic",
             multiview = TRUE,
             foo = playRectInput(playState, prompt = prompt,
                                            scales = scales))
{
    force(foo)
    if (is.null(foo)) return(NULL)
    if (is.null(foo$coords)) return(NULL)
    coords <- foo$coords
    data <- xyData(playState, space = foo$space)
    if (length(data$x) == 0) {
        gmessage.error(paste("Sorry, can not guess the data point coordinates.",
                             "Please contact the maintainer with suggestions."))
        return(NULL)
    }
    ## convert to numeric (same as xyCoords())
    numdata <- data
    ## assume if not numeric, then not a matrix
    if (!is.numeric(data$x))
        numdata$x <- as.numeric(data$x)
    if (!is.numeric(data$y))
        numdata$y <- as.numeric(data$y)
    ## convert to log scale if necessary
    numdata <- dataCoordsToSpaceCoords(playState, numdata)

    which <- NULL
    pos <- NULL
    if (foo$is.click) {
        x <- coords$x[1]
        y <- coords$y[1]
        ppxy <- playDo(playState,
                       quote(list(
                            lx=convertX(unit(x, "native"), "points", TRUE),
                            ly=convertY(unit(y, "native"), "points", TRUE),
                            px=convertX(unit(numdata$x, "native"), "points", TRUE),
                            py=convertY(unit(numdata$y, "native"), "points", TRUE))),
                       space=foo$space)
        pdists <- with(ppxy, sqrt((px - lx)^2 + (py - ly)^2))
        if (min(pdists, na.rm = TRUE) > 18)
            which <- integer(0)
        else {
            which <- which.min(pdists)
            pos <- with(ppxy, getTextPosition(x = lx - px[which],
                                              y = ly - py[which]))
        }
    }
    else {
        ## drag
        ok <- TRUE
        if (!foo$yOnly)
            ok <- (min(coords$x) <= numdata$x) & (numdata$x <= max(coords$x))
        if (!foo$xOnly)
            ok <- ok & (min(coords$y) <= numdata$y) & (numdata$y <= max(coords$y))
        which <- which(ok)
    }
    ## locations of actual points clicked
    x <- data$x[which]
    y <- data$y[which]
    ## account for multiple points (matrix values of data$x, data$y)
    n <- min(NROW(data$x), NROW(data$y))
    which <- unique(which %% n)
    which[which == 0] <- n
    if (multiview) {
        x <- if (is.matrix(data$x)) data$x[which,,drop=FALSE] else data$x[which]
        y <- if (is.matrix(data$y)) data$y[which,,drop=FALSE] else data$y[which]
    }
    subscripts <- data$subscripts[which]
    if (is.null(subscripts)) subscripts <- which
    c(list(subscripts = subscripts, which = which,
           x = x, y = y,
           pos = pos, is.click = foo$is.click),
      foo)
}

playGetIDs <-
    function(playState = playDevCur(),
             type = c("labelled", "brushed"),
             labels = FALSE)
{
    type <- match.arg(type, several.ok = TRUE)
    ids.brushed <- unlist(playState$linked$ids)
    ids.labelled <- do.call(rbind, playState$ids)$subscripts
    ids <- NULL
    if ("labelled" %in% type) ids <- ids.labelled
    if ("brushed" %in% type) ids <- c(ids, ids.brushed)
    if (length(ids) > 0)
        ids <- unique(sort(ids))
    if (labels) playState$labels[ids] else ids
}

playSetIDs <-
    function(playState = playDevCur(),
             value,
             type = "brushed",
             space = "plot",
             add = FALSE,
             redraw = NA,
             pos = 1)
{
    playStoreUndo(playState)
    type <- match.arg(type, c("labelled", "brushed"))
    if (is.logical(value))
        value <- which(value)
    if (type == "brushed") {
        if (add == FALSE) { ## replace
            if (length(playState$linked$ids)) ## exists
                if (is.na(redraw)) redraw <- TRUE
            playState$linked$ids <- list()
        }
        i <- length(playState$linked$ids) + 1
        playState$linked$ids[[i]] <- value
        if (is.na(redraw)) {
            ## draw without a full redraw
            drawLinkedLocal(playState)
            updateLinkedSubscribers(playState)
        }
    }
    if (type == "labelled") {
        if (add == FALSE) { ## replace
            if (length(playState$ids)) ## exists
                if (is.na(redraw)) redraw <- TRUE
            playState$ids <- list()
        }
        ids.new <- data.frame(subscripts = value, pos = pos)
        i <- length(playState$ids) + 1
        playState$ids[[i]] <- ids.new
        names(playState$ids)[i] <- space
        if (is.na(redraw)) {
            ## draw without a full redraw
            drawLabelsInSpace(playState, subscripts = value,
                              space = space, pos = pos)
            #drawLabels(playState)
        }
    }
    if (isTRUE(redraw)) {
        ## full redraw
        playReplot(playState)
        ## redraw linked plots
        if (type == "brushed")
            updateLinkedSubscribers(playState, redraw = TRUE)
    }
    ## update other tool states
    updateAnnotationActionStates(playState)
    updateIdentifyActionStates(playState)
    invisible()
}

playClear <-
    function(playState = playDevCur(),
             type = c("annotations", "labelled", "brushed"),
             redraw = TRUE)
{
    playStoreUndo(playState)
    type <- match.arg(type, several.ok = TRUE)
    ## remove types that are empty
    type <- c(if (("labelled" %in% type) &&
                  length(playState$ids)) "labelled",
              if (("annotations" %in% type) &&
                  length(playState$annotations)) "annotations",
              if (("brushed" %in% type) &&
                  length(playState$linked$ids)) "brushed"
              )
    if (length(type) == 0) return()
    for (x in type) {
        if (x == "labelled") {
            playState$ids <- list()
        } else if (x == "annotations") {
            playState$annotations <- list()
        } else if (x == "brushed") {
            playState$linked$ids <- list()
        }
    }
    ## redraw
    if (redraw) {
        playReplot(playState)
        ## update linked plots
        if ("brushed" %in% type)
            updateLinkedSubscribers(playState, redraw = TRUE)
    }
    updateAnnotationActionStates(playState)
    updateIdentifyActionStates(playState)
    invisible()
}

playUndo <- function(playState = playDevCur())
{
    if (isBasicDeviceMode(playState)) {
        ## basic device mode: only one stored display
        redoPlot <- recordPlot()
        try(replayPlot(playState$tmp$recorded.plot))
        generateSpaces(playState)
        playState$tmp$recorded.plot <- redoPlot
        return(invisible())
    }
    i <- length(playState$tmp$undoStack)
    if (i == 0) return()
    state <- playState$tmp$undoStack[[i]]
    anyLinked <- !identical(state$brushed, playState$linked$ids)
    length(playState$tmp$undoStack) <- (i - 1)
    playState$ids <- state$ids
    playState$annotations <- state$annotations
    playState$linked$ids <- state$brushed
    ## redraw
    playReplot(playState)
    if (anyLinked)
        updateLinkedSubscribers(playState, redraw = TRUE)
    invisible()
}

playStoreUndo <- function(playState = playDevCur())
{
    if (isBasicDeviceMode(playState)) {
        ## basic device mode: only one stored display
        playState$tmp$recorded.plot <- try(recordPlot())
        return(invisible())
    }
    state <- list()
    state$ids <- playState$ids
    state$annotations <- playState$annotations
    state$brushed <- playState$linked$ids
    i <- length(playState$tmp$undoStack) + 1
    playState$tmp$undoStack[[i]] <- state
    if (i > playwith.getOption("undo.levels"))
        playState$tmp$undoStack <- playState$tmp$undoStack[-1]
    updateAnnotationActionStates(playState)
    invisible()
}

drawLabels <- function(playState, return.code = FALSE)
{
    playDevSet(playState)
    theCode <- expression()
    ## group by space
    spaces <- names(playState$ids)
    for (space in unique(spaces)) {
        items <- playState$ids[spaces == space]
        idInfo <- do.call(rbind, items)
        expr <- drawLabelsInSpace(playState, subscripts = idInfo$subscripts,
                           space = space, pos = idInfo$pos,
                           return.code = return.code)
        if (return.code)
            theCode <- c(theCode, expr)
    }
    theCode
}

drawLabelsInSpace <- function(playState, subscripts, space = "plot",
                              pos = 1, return.code = FALSE)
{
    data <- xyCoords(playState, space=space)
    if (length(data$x) == 0) return()
    if (length(data$y) == 0) return()
    ## convert to log scale if necessary
    data <- dataCoordsToSpaceCoords(playState, data)
    if (!is.null(data$subscripts)) {
        ## 'data' is a subset given by data$subscripts,
        ## so need to find which ones match the label subscripts
        which <- match(subscripts, data$subscripts, 0)
    } else {
        ## 'data' (x and y) is the whole dataset
        which <- subscripts
    }
    x <- if (is.matrix(data$x)) data$x[which,] else data$x[which]
    y <- if (is.matrix(data$y)) data$y[which,] else data$y[which]
    nvar <- max(NCOL(data$x), NCOL(data$y))
    labels <- playState$labels[subscripts]
    pos <- rep(pos, length = length(labels))
    ## if each data point has multiple locations, replicate labels
    labels <- rep(labels, each = nvar)
    pos <- rep(pos, each = nvar)
    offset <- as.numeric(playState$label.offset)
    annots <- expression()
    for (i in seq_along(labels)) {
        annots[[i]] <- call("panel.usertext", x[i], y[i],
                            labels[i], pos = pos[i])
        if (offset != 0.5)
            annots[[i]]$offset <- offset
    }
    playDo(playState, annots, space = space,
           return.code = return.code)
}

drawLinkedLocal <- function(playState, return.code = FALSE)
{
    ## draw linked brushed points
    theCode <- expression()
    subscripts <- unlist(playState$linked$ids)
    if (length(subscripts) == 0) return(theCode)
    subscripts <- unique(sort(subscripts))
    playDevSet(playState)
    for (space in playState$spaces) {
        data <- xyCoords(playState, space = space)
        if (length(data$x) == 0) next
        if (length(data$y) == 0) next
        ## convert to log scale if necessary
        data <- dataCoordsToSpaceCoords(playState, data)
        if (!is.null(data$subscripts)) {
            ## 'data' is a subset given by data$subscripts,
            ## so need to find which ones match the label subscripts
            which <- match(subscripts, data$subscripts, 0)
        } else {
            ## 'data' (x and y) is the whole dataset
            which <- subscripts
        }
        x <- if (is.matrix(data$x)) data$x[which,,drop=FALSE] else data$x[which]
        y <- if (is.matrix(data$y)) data$y[which,,drop=FALSE] else data$y[which]
        if (length(x) == 0) next
        annot <- call("panel.brushpoints", x, y)
        ## special case: parallel -- draw lines not points
        if (playState$callName == "parallel") {
            ## use NAs to achieve breaks in lines
            x <- t(cbind(x, NA))
            y <- t(cbind(y, NA))
            annot <- call("panel.brushlines", x, y)
        }
        expr <- playDo(playState, annot, space = space,
                       return.code = return.code)
        if (return.code)
            theCode <- c(theCode, expr)
    }
    theCode
}

updateLinkedSubscribers <- function(playState, redraw = FALSE)
{
    whichDead <- NULL
    for (i in seq_along(playState$linked$subscribers)) {
        otherPlayState <- playState$linked$subscribers[[i]]
        if (!identical(otherPlayState$ID, playState$ID)) {
            ## first check that this subscriber is still alive
            if (!inherits(otherPlayState$win, "GtkWindow")) {
                whichDead <- c(whichDead, i)
                next
            }
            ## trigger draw / redraw
            if (redraw) playReplot(otherPlayState)
            else drawLinkedLocal(otherPlayState)
        }
    }
    if (length(whichDead))
        playState$linked$subscribers <-
            playState$linked$subscribers[-whichDead]
}

playUnlink <- function(playState = playDevCur())
{
    oldlinked <- playState$linked
    newlinked <- new.env(parent = baseenv())
    newlinked$ids <- playState$linked$ids
    newlinked$subscribers <- list(playState)
    playState$linked <- newlinked
    ## remove self from (old) list of subscribers
    for (i in seq_along(oldlinked$subscribers)) {
        otherPlayState <- oldlinked$subscribers[[i]]
        if (identical(otherPlayState$ID, playState$ID)) {
            oldlinked$subscribers <- oldlinked$subscribers[-i]
            break
        }
    }
    ## update action states
    for (x in oldlinked$subscribers)
        updateGlobalActions(x)
    updateGlobalActions(playState)
}



## COPIED FROM LATTICE
getTextPosition <- function(x, y)
    ## returns position 1: below, 2: left, 3: above, 4: right (w.r.t
    ## origin).  Used as a tool in panel.identify.
{
    a <- abs(c(x, y))
    if (y <= 0 && a[1] <= -y) 1
    else if (x <= 0 && a[2] <= -x) 2
    else if (y >= 0 && a[1] <= y) 3
    else if (x >= 0 && a[2] <= x) 4
}
