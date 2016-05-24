## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer


updateAnnotationActions <- function(playState)
{
    drawAnnotations(playState)
    updateAnnotationActionStates(playState)
}

updateAnnotationActionStates <- function(playState)
{
    aGroup <- playState$actionGroups[["PlotActions"]]
    ## UndoAnnotation
    showUndo <- (length(playState$tmp$undoStack) > 0)
    if (isBasicDeviceMode(playState))
        showUndo <- !is.null(playState$tmp$recorded.plot)
    aGroup$getAction("UndoAnnotation")$setSensitive(showUndo)
    ## Clear
    ## in basic device mode do not know the call; it cannot be redrawn
    showClear <- !isBasicDeviceMode(playState)
    showClear <- showClear && ((length(playState$ids) > 0) ||
                               (length(playState$annotations) > 0) ||
                               (length(playState$linked$ids) > 0))
    aGroup$getAction("Clear")$setSensitive(showClear)
    ## EditAnnotations
    showEdit <- !isBasicDeviceMode(playState)
    showEdit <- showEdit && (length(playState$annotations) > 0)
    aGroup$getAction("EditAnnotations")$setSensitive(showEdit)
}

clear_handler <- function(widget, playState)
{
    types <- c(
               if (length(playState$ids) > 0) "labelled",
               if (length(playState$linked$ids) > 0) "brushed",
               if (length(playState$annotations) > 0) "annotations"
               )
    if (length(types) == 0) { widget$hide(); return() }
    clear.types <- types
    if (length(types) > 1) {
        clear.types <- NULL
        widg <- gcheckboxgroup(types, checked=TRUE)
        result <- gbasicdialog(title="Clear what?", widget=widg,
                               handler=function(...)
                               clear.types <<- svalue(widg) )
        playState$win$present()
        if (!isTRUE(result)) return()
    }
    playClear(playState, type = clear.types)
}

undo.annotation_handler <- function(widget, playState)
    playUndo(playState)

edit.annotations_handler <- function(widget, playState)
{
    annotTxt <- deparse(playState$annotations, width.cutoff = 42,
                        control = playwith.getOption("deparse.options"))
    annotTxt <- paste(annotTxt, collapse = "\n")
    ## TODO: deparse / parse in a more readable form
    repeat {
        newTxt <-
            guiTextInput(annotTxt, title="Edit annotations",
                         prompt=paste("Make sure you keep this structure!",
                         "(a list of expressions, named by space)", sep="\n"),
                         accepts.tab=FALSE)
        if (is.null(newTxt)) break
        annotTxt <- newTxt
        tmp <- tryCatch(eval(parse(text=annotTxt)), error=function(e)e)
        ## check whether there was a syntax error
        if (inherits(tmp, "error")) {
            gmessage.error(conditionMessage(tmp))
        } else {
            playStoreUndo(playState)
            playState$annotations <- tmp
            playReplot(playState)
            ## update other tool states
            updateAnnotationActionStates(playState)
            break
        }
    }
    playState$win$present()
}

rectCore <- function(playState, foo)
{
    pageAnnotation <- isTRUE(playState$page.annotation)
    if (is.null(foo)) return()
    if (is.null(foo$coords)) pageAnnotation <- TRUE
    if (foo$is.click) return()
    space <- foo$space
    if (pageAnnotation) space <- "page"
    myXY <- if (space == "page") foo$ndc else foo$coords
    myXY$x <- signif(myXY$x, 7)
    myXY$y <- signif(myXY$y, 7)
    annot <- with(myXY, call("panel.rect",
                             min(x), min(y), max(x), max(y)))
    ## draw it and store it
    playAnnotate(playState, annot, space = space)
}

lineCore <- function(playState, foo)
{
    pageAnnotation <- isTRUE(playState$page.annotation)
    if (is.null(foo)) return()
    if (is.null(foo$coords)) pageAnnotation <- TRUE
    if (foo$is.click) return()
    space <- foo$space
    if (pageAnnotation) space <- "page"
    myXY <- if (space == "page") foo$ndc else foo$coords
    myXY$x <- signif(myXY$x, 7)
    myXY$y <- signif(myXY$y, 7)
    annot <- with(myXY, call("panel.segments",
                             x[1], y[1], x[2], y[2]))
    ## draw it and store it
    playAnnotate(playState, annot, space = space)
}

arrowCore <- function(playState, foo)
{
    pageAnnotation <- isTRUE(playState$page.annotation)
    if (is.null(foo)) return()
    if (is.null(foo$coords)) pageAnnotation <- TRUE
    if (foo$is.click) return()
    space <- foo$space
    if (pageAnnotation) space <- "page"
    myXY <- if (space == "page") foo$ndc else foo$coords
    myXY$x <- signif(myXY$x, 7)
    myXY$y <- signif(myXY$y, 7)
    annot <- with(myXY, call("panel.arrows",
                             x[1], y[1], x[2], y[2]))
    arrow <- playState$arrow
    ## each of these may be NULL
    annot$angle <- arrow$angle
    annot$length <- arrow$length
    annot$unit <- arrow$unit
    annot$type <- arrow$type
    annot$ends <- arrow$ends
    annot$code <- arrow$code
    ## draw it and store it
    playAnnotate(playState, annot, space = space)
}

annotateCore <- function(playState, foo)
{
    pageAnnotation <- isTRUE(playState$page.annotation)
    if (is.null(foo)) return()
    if (is.null(foo$coords)) pageAnnotation <- TRUE
    space <- foo$space
    if (pageAnnotation) space <- "page"
    absXY <- foo$ndc
    myXY <- if (space == "page") foo$ndc else foo$coords
    myXY$x <- signif(myXY$x, 7)
    myXY$y <- signif(myXY$y, 7)
    if (foo$is.click) {
        myX <- myXY$x[1]
        myY <- myXY$y[1]
    }

    playFreezeGUI(playState)
    on.exit(playThawGUI(playState))

    ## pop up dialog to create label
    dialog <- gwindow(title = "New annotation")
    wingroup <- ggroup(horizontal = FALSE, container = dialog)
    wid <- list()

    ## TEXT AND JUST
    labgroup <- gframe("Label text and position",
                       horizontal = FALSE, container = wingroup)
    wid$label <- gtext(width = 200, height = 50, container = labgroup)
    wid$label.expr <- gcheckbox("plotmath", container = labgroup)
    just.hgroup <- ggroup(horizontal = TRUE, container = labgroup)
    glabel(if (foo$is.click) "Position of text \n relative to click:" else
           "Justification of \n text inside box:", container = just.hgroup)
    ## this grid of checkboxes is conceptually a radio button group
    justw <- list()
    just_handler <- function(h, ...) {
        if (svalue(h$obj) == FALSE) return()
        ## turn all the others off
        for (j in names(justw)) {
            if (j != h$action)
                svalue(justw[[j]]) <- FALSE
        }
    }
    lay <- glayout(container=just.hgroup, spacing=2)
    lay[1,1] <- justw$lt <- gcheckbox("topleft",  handler = just_handler, action = "lt")
    lay[1,2] <- justw$ct <- gcheckbox("top    ",  handler = just_handler, action = "ct")
    lay[1,3] <- justw$rt <- gcheckbox("topright", handler = just_handler, action = "rt")
    lay[2,1] <- justw$lc <- gcheckbox("left   ",  handler = just_handler, action = "lc")
    lay[2,2] <- justw$cc <- gcheckbox("centre ",  handler = just_handler, action = "cc")
    lay[2,3] <- justw$rc <- gcheckbox("right  ",  handler = just_handler, action = "rc")
    lay[3,1] <- justw$lb <- gcheckbox("botleft",  handler = just_handler, action = "lb")
    lay[3,2] <- justw$cb <- gcheckbox("bottom ",  handler = just_handler, action = "cb")
    lay[3,3] <- justw$rb <- gcheckbox("botright", handler = just_handler, action = "rb")
    svalue(justw$cc) <- TRUE
    visible(lay) <- TRUE

    offsetgroup <- ggroup(horizontal = TRUE, container = labgroup)
    glabel(if (foo$is.click) "Offset from point (chars):"
    else "Offset from box edge", container = offsetgroup)
    wid$offset <- gedit(toString(playState$label.offset), width = 4,
                        coerce.with = as.numeric, container = offsetgroup)
    if (!foo$is.click) {
        ## it was a drag, defining a rectangle
        ## option to draw box border
        wid$drawbox <- gcheckbox("Draw box", container = offsetgroup)
        ## TODO: fit to box?
    }
    focus(wid$label) <- TRUE

    ## STYLE
    user.text <- trellis.par.get("user.text")
    if (is.null(eval(user.text)))
        user.text <- trellis.par.get("add.text")
    stylegroup <- gframe("Style",
                         horizontal=FALSE, container=wingroup)
    lay <- glayout(container=stylegroup, spacing=2)
    ## style settings widgets
    wid$col <- gdroplist(palette(), selected = 0, editable = TRUE)
    size(wid$col) <- c(100, -1)
    wid$cex <- gedit("1.0", width = 4, coerce.with = as.numeric)
    wid$lineheight <- gedit("1.0", width = 4, coerce.with = as.numeric)
    wid$rot <- gdroplist(c("-90","-45","-30","0","30","45","90"), selected = 4,
                         editable = TRUE, coerce.with = as.numeric)
    size(wid$rot) <- c(60, -1)
    if (!is.null(user.text$col)) svalue(wid$col) <- user.text$col
    if (!is.null(user.text$cex)) svalue(wid$cex) <- user.text$cex
    if (!is.null(user.text$lineheight))
        svalue(wid$lineheight) <- user.text$lineheight
    lay[1,1] <- "Text color:"
    lay[1,2] <- wid$col
    lay[1,3] <- " Scale:"
    lay[1,4] <- wid$cex
    lay[2,1] <- "Lineheight:"
    lay[2,2] <- wid$lineheight
    lay[2,3] <- " Rotation:"
    lay[2,4] <- wid$rot
    visible(lay) <- TRUE
    wid$set.defaults <- gcheckbox("Set as default label style",
                                  container = stylegroup)
    glabel("For more control, try Style Settings from the Style menu.",
           container = stylegroup)

    ## store plot display for fast redrawing
    originalPlot <- try(recordPlot())
    showingPreview <- FALSE

    annot_handler <- function(h, ...)
    {
        ## note: playState is accessed from the function environment

        ## TEXT AND JUST
        argExpr <- function(wid, expr.wid) {
            newVal <- svalue(wid)
            if (newVal == "") return(NULL)
            if (svalue(expr.wid))
                newVal <- try(parse(text=newVal, srcfile=NULL))
            if (inherits(newVal, "try-error")) return(NULL)
            newVal
        }
        isExpr <- svalue(wid$label.expr)
        labelVal <- argExpr(wid$label, wid$label.expr)
        if (!isExpr) {
            labelVal <- gsub("\\", "\\\\", labelVal, fixed=TRUE)
        }
        ## find which checkbox was selected
        just <- c("c", "c")
        for (j in names(justw))
            if (svalue(justw[[j]])) just <- strsplit(j, NULL)[[1]]
        just[1] <- switch(just[1], l="left", r="right", c="centre")
        just[2] <- switch(just[2], t="top", b="bottom", c="centre")

        ## COORDINATES
        if (foo$is.click == FALSE) {
            ## it was a drag, defining a rectangle
            ## choose side of rect to align to
            x.leftright <- myXY$x[order(absXY$x)]
            y.bottop <- myXY$y[rev(order(absXY$y))]
            myX <- switch(just[1],
                          left = x.leftright[1],
                          right = x.leftright[2],
                          mean(x.leftright))
            myY <- switch(just[2],
                          bottom = y.bottop[1],
                          top = y.bottop[2],
                          mean(y.bottop))
            ## justification flipped (inside rect, not around point)
            just[1] <- switch(just[1],
                              left="right", right="left", centre="centre")
            just[2] <- switch(just[2],
                              top="bottom", bottom="top", centre="centre")
        }
        ## have to add offset manually if text is placed on a corner
        if ((svalue(wid$offset) != 0) && all(just != "centre")) {
            pad <- svalue(wid$offset)
            pad <- playDo(playState,
                          convertWidth(unit(pad, "char"), "native", TRUE),
                          space = space)
            if (just[1] == "left") pad <- 0 - pad
            myX <- myX + pad
            pad <- svalue(wid$offset)
            pad <- playDo(playState,
                          convertHeight(unit(pad, "char"), "native", TRUE),
                          space = space)
            if (just[2] == "bottom") pad <- 0 - pad
            myY <- myY + pad
        }

        myX <- signif(myX, 7)
        myY <- signif(myY, 7)

        ## CREATE THE CALL
        annot <- call("panel.usertext", myX, myY, labelVal)
        if (!all(just == "centre")) {
            if (any(just == "centre")) {
                if ("bottom" %in% just) annot$pos <- 1
                if ("left" %in% just) annot$pos <- 2
                if ("top" %in% just) annot$pos <- 3
                if ("right" %in% just) annot$pos <- 4
                if (svalue(wid$offset) != 0.5)
                    annot$offset <- svalue(wid$offset)
            } else {
                adj <- c(0.5, 0.5)
                if ("bottom" %in% just) adj[2] <- 1
                if ("left" %in% just) adj[1] <- 1
                if ("top" %in% just) adj[2] <- 0
                if ("right" %in% just) adj[1] <- 0
                annot$adj <- adj
            }
        }
        if (svalue(wid$set.defaults)) {
            user.text$col <- svalue(wid$col)
            user.text$cex <- svalue(wid$cex)
            user.text$lineheight <- svalue(wid$lineheight)
            trellis.par.set(user.text = user.text)
        }
        if (!identical(svalue(wid$col), user.text$col))
            annot$col <- svalue(wid$col)
        if (!identical(svalue(wid$cex), user.text$cex))
            annot$cex <- svalue(wid$cex)
        if (!identical(svalue(wid$lineheight), user.text$lineheight))
            annot$lineheight <- svalue(wid$lineheight)
        if (!identical(svalue(wid$rot), 0))
            annot$srt <- svalue(wid$rot)

        ## draw box?
        if (!foo$is.click && svalue(wid$drawbox)) {
            dobox <- call("panel.rect", x = mean(myXY$x), y = mean(myXY$y),
                          width = abs(diff(myXY$x)), height = abs(diff(myXY$y)))
            annot <- call("{", dobox, annot)
        }

        ## need to redraw all annotations if style changed
        if (svalue(wid$set.defaults) &&
            (length(playState$annotations) > 0) &&
            !isBasicDeviceMode(playState))
        {
            playReplot(playState)
            originalPlot <- try(recordPlot())
            showingPreview <<- FALSE
        }

        if (showingPreview) {
            ## remove preview (i.e. redraw original plot)
            result <- try(replayPlot(originalPlot))
            ## may fail if engine.display.list is off
            if (inherits(result, "try-error")) {
                playReplot(playState)
            } else {
                generateSpaces(playState)
            }
        }

        if (h$action == "preview") {
            ## echo annotation code to console
            message(paste(deparse(annot), collapse="\n"))
            ## draw it
            playDo(playState, annot, space = space)
            showingPreview <<- TRUE
            return()
        }
        ## draw it and store it
        playAnnotate(playState, annot, space = space)

        dispose(h$obj)
        playState$win$present()
    }

    buttgroup <- ggroup(container=wingroup)
    addSpring(buttgroup)
    okbutt <- gbutton("OK", handler=annot_handler,
                      action="ok", container=buttgroup)
    prebutt <- gbutton("Preview", handler=annot_handler,
                       action="preview", container=buttgroup)
    canbutt <- gbutton("Cancel", handler=function(h, ...) {
        if (showingPreview) {
            ## remove preview (i.e. redraw original plot)
            result <- try(replayPlot(originalPlot))
            if (inherits(result, "try-error")) {
                playReplot(playState)
            } else {
                generateSpaces(playState)
            }
        }
        dispose(h$obj)
    }, container=buttgroup)
    size(okbutt) <- size(prebutt) <- size(canbutt) <- c(80, 30)
    #defaultWidget(okbutt) <- TRUE
}

legend_handler <- function(widget, playState)
{
    ## TODO place legend
}

