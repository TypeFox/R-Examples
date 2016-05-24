## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

initClickActions <- function(playState)
{
    ## add click handler to device
    if (is.null(playState$widgets$buttonPressSignal)) {
        playState$widgets$buttonPressSignal <-
            gSignalConnect(playState$widgets$drawingArea,
                           "button-press-event",
                           device.click_handler, data=playState)
    }
    if (!is.null(playState$click.mode)) {
        ## set initial click.mode
        vals <- clickModeValues()
        val <- vals[[playState$click.mode]]
        if (!is.null(val)) {
            aGroup <- playState$actionGroups[["PlotActions"]]
            try(aGroup$getAction("Zoom")$setCurrentValue(val),
                silent = TRUE)
        }
        playState$click.mode <- NULL
    }
}

clickmode.change_handler <- function(action, current, playState)
{
    if (current["active"]) {
        playState$tmp$click.mode <- gtkActionGetName(current)
        try(updateClickActions(playState))
    }
    return()
}

updateClickActions <- function(playState)
{
    if (is.null(playState$tmp$click.mode))
        playState$tmp$click.mode <- "Zoom"
    curs <- switch(playState$tmp$click.mode,
                   Zoom = "crosshair",
                   Pan = "fleur",
                   Identify = "hand1",
                   Brush = "circle",
                   Annotation = "xterm",
                   "left_ptr") ## (otherwise)
    cursor <- gdkCursorNew(GdkCursorType[curs])
    playState$widgets$drawingArea$getWindow()$setCursor(cursor)
    ## work out which actions are possible on current plot
    hasArgs <- playState$accepts.arguments
    isLatt <- playState$is.lattice
    isSplom <- (playState$callName %in% c("splom"))
    isLatt3D <- isLatt && !is.null(playState$trellis$panel.args.common$scales.3d)
    hasPanels <- isLatt && (length(playState$tmp$currentLayout) > 1)
    isBase <- !isLatt && is.null(playState$viewport)
    isBaseMulti <- isBase && any(par("mfrow") > 1)
    canIdent <- isTRUE(playState$tmp$identify.ok)
    actions <- list()
    actions$nav2D <- (hasArgs && !isLatt3D && !isSplom)
    actions$nav3D <- (isLatt3D)
    actions$ident <- (canIdent)
    playState$tmp$ok.actions <- actions
    ## set default statusbar message
    modeOK <- playState$tmp$click.mode
    if ((modeOK == "Pan") && actions$nav3D)
        modeOK <- "Rotate"
    if ((modeOK %in% c("Zoom", "Pan")) &&
        (actions$nav2D == FALSE))
        modeOK <- "Coords"
    if ((modeOK %in% c("Identify", "Brush")) &&
        (actions$ident == FALSE))
        modeOK <- "Coords"
    msg <- switch(modeOK,
                  Zoom = paste("Drag to zoom (hold Shift to constrain)",
                  "Click for coordinates", sep = ", "),
                  Pan = paste("Drag to scroll (hold Shift to constrain)",
                  "Click for coordinates", sep = ", "),
                  Rotate = "Drag to rotate (hold Shift to constrain)",
                  Identify = "Click or drag to identify points",
                  Brush = paste("Click or drag to brush points",
                  "(hold Shift to add, not replace)"),
                  Annotation = "Click or drag to place text",
                  Arrow = "Drag to draw an arrow (hold Shift to constrain)",
                  Line = "Drag to draw a line (hold Shift to constrain)",
                  Rect = "Drag to draw a rectangle",
                  Coords = "Click for coordinates") ## fallback
    ## Zoom actions are always accessible, if possible:
    if (actions$nav2D || actions$nav3D) {
        if (modeOK != "Zoom")
            msg <- paste(msg, "Ctrl-drag to zoom", sep = ", ")
        msg <- paste(msg, "Ctrl-click to zoom out", sep = ", ")
    }
    msg <- paste(msg, "Right-click for more", sep = ", ")
    playState$widgets$statusbar$pop(1)
    playState$widgets$statusbar$push(1, msg)
}

device.click_handler <- function(widget, event, playState)
{
    if (!isTRUE(playState$tmp$plot.ready)) return(FALSE)
    ## bail out if another tool is handling the click (this ok?)
    if (isTRUE(playState$tmp$now.interacting)) return(FALSE)
    playDevSet(playState)
    if (playState$tmp$need.reconfig) generateSpaces(playState)
    x <- event$x
    y <- event$y
    ## work out which actions are relevant to the plot
    actions <- playState$tmp$ok.actions
    modeOK <- playState$tmp$click.mode
    if ((modeOK == "Pan") && actions$nav3D)
        modeOK <- "Rotate"
    if ((modeOK %in% c("Zoom", "Pan")) &&
        (actions$nav2D == FALSE))
        modeOK <- "Coords"
    if ((modeOK %in% c("Identify", "Brush")) &&
        (actions$ident == FALSE))
        modeOK <- "Coords"
    pageOK <- (modeOK %in% c("Annotation", "Arrow", "Line", "Rect"))
    isCtrlClick <- (as.flag(event$state) & GdkModifierType["control-mask"])
    if (.Platform$OS.type == "unix") {
        isAltClick <- (as.flag(event$state) & GdkModifierType["mod1-mask"])
    } else {
        isAltClick <- ((as.flag(event$state) & GdkModifierType["mod1-mask"]) ||
                       (as.flag(event$state) & GdkModifierType["mod2-mask"]))
    }
    isShiftClick <- (as.flag(event$state) & GdkModifierType["shift-mask"])
    isPlainClick <- !isCtrlClick && !isAltClick && !isShiftClick
    ## take actions
    ## Alt-click is treated as a right-click, intended for MacOS (untested)
    if ((event$button == 1) && !isAltClick) {
        ## standard (left) mouse button
        dragShape <- "rect"
        if (!isAltClick) {
            if (modeOK %in% c("Pan", "Rotate", "Arrow", "Line"))
                dragShape <- "line"
        }
        scales <- "dynamic"
        if (modeOK %in% c("Identify", "Brush"))
            scales <- c("x", "y")
        ## handle click or drag
        foo <- playClickOrDrag(playState, x0=x, y0=y,
                               shape=dragShape, scales = scales)
        if (is.null(foo)) {
            ## drag went off device
            coordsCore(playState, NULL)
            return(FALSE)
        }
        if (is.null(foo$coords) && !pageOK) {
            ## click/drag outside of a defined space
            coordsCore(playState, NULL)
            return(FALSE)
        }
        if (foo$is.click) {
            coordsCore(playState, foo)
        } else {
            coordsCore(playState, NULL)
        }
        ## standard Ctrl-click actions
        if (isCtrlClick) {
            if (actions$nav2D) {
                ## 2D Zoom
                if (foo$is.click) {
                    zoomoutCore(playState, foo)
                } else {
                    ## drag
                    zoomCore(playState, foo)
                }
            }
            if (actions$nav3D) {
                ## 3D Zoom
                if (foo$is.click) {
                    zoomout3DCore(playState, foo)
                } else {
                    ## drag
                    zoom3DCore(playState, foo)
                }
            }
        } else {
            ## plain click: normal actions
            if (modeOK == "Zoom") {
                if (!foo$is.click)
                    zoomCore(playState, foo)
            }
            if (modeOK == "Pan") {
                panCore(playState, foo)
            }
            if (modeOK == "Rotate") {
                if (!foo$is.click)
                    rotate3DCore(playState, foo)
            }
            if (modeOK == "Identify") {
                identifyCore(playState, foo)
            }
            if (modeOK == "Brush") {
                brushCore(playState, foo, add = isShiftClick)
            }
            if (modeOK == "Annotation") {
                annotateCore(playState, foo)
            }
            if (modeOK == "Arrow") {
                arrowCore(playState, foo)
            }
            if (modeOK == "Line") {
                lineCore(playState, foo)
            }
            if (modeOK == "Rect") {
                rectCore(playState, foo)
            }
        }
    } else {
        coordsCore(playState, NULL)
    }
    if (event$button == 2) {
        ## middle mouse button click: zoom to fit
        if (actions$nav2D || actions$nav3D)
            zoomfit_handler(NULL, playState)
    }
    if ((event$button == 3) || isAltClick) {
        ## right mouse button or alt-click
        foo <- playClickOrDrag(playState, x0=x, y0=y,
                               shape="rect")
        if (is.null(foo))
            return(FALSE)
        ## pop up context menu
        contextCore(playState, foo, event = event)
    }
    return(FALSE)
}

coordsCore <- function(playState, foo) {
    coords <- foo$coords
    ## convert from log scale if necessary
    coords <- spaceCoordsToDataCoords(playState, coords)
    if (!is.null(coords)) {
        x <- format(coords$x, digits=3)
        y <- format(coords$y, digits=3)
        coordsTxt <- paste("<tt>x:", x, ", y:", y, " </tt>", sep="")
        playState$widgets$coordsLabel$setMarkup(coordsTxt)
        playState$widgets$coordsLabel["visible"] <- TRUE
    } else {
        playState$widgets$coordsLabel["visible"] <- FALSE
    }
}

panCore <- function(playState, foo)
{
    xdiff <- diff(foo$coords$x)
    ydiff <- diff(foo$coords$y)
    ## update axis scales
    if (!foo$yOnly) {
        xlim <- rawXLim(playState, space=foo$space)
        rawXLim(playState) <- xlim - xdiff
    }
    if (!foo$xOnly) {
        ylim <- rawYLim(playState, space=foo$space)
        rawYLim(playState) <- ylim - ydiff
    }
    playReplot(playState)
}

zoomCore <- function(playState, foo)
{
    xlim <- range(foo$coords$x)
    ylim <- range(foo$coords$y)
    ## update axis scales, reverse if needed
    if (!foo$yOnly) {
        if (is.unsorted(rawXLim(playState, space=foo$space)))
            xlim <- rev(xlim)
        rawXLim(playState) <- xlim
    }
    if (!foo$xOnly) {
        if (is.unsorted(rawYLim(playState, space=foo$space)))
            ylim <- rev(ylim)
        rawYLim(playState) <- ylim
    }
    playReplot(playState)
}

zoomoutCore <- function(playState, foo)
{
    nav.x <- !isTRUE(foo$yOnly)
    #nav.y <- !isTRUE(foo$xOnly)
    ## in time.mode, only zoom along x-axis
    nav.y <- (!isTRUE(playState$time.mode) &&
              is.null(playState$time.vector))
    ## find existing scales
    xlim <- rawXLim(playState, space=foo$space)
    ylim <- rawYLim(playState, space=foo$space)
    ## centre on click location
    xlim <- (xlim - mean(xlim)) + mean(foo$coords$x)
    ylim <- (ylim - mean(ylim)) + mean(foo$coords$y)
    ## zoom out: make range twice the size
    if (nav.x) xlim <- xlim + diff(xlim) * c(-0.5, 0.5)
    if (nav.y) ylim <- ylim + diff(ylim) * c(-0.5, 0.5)
    ## this converts from raw numeric to original format (including unlog)
    if (nav.x) rawXLim(playState) <- xlim
    if (nav.y) rawYLim(playState) <- ylim
    playReplot(playState)
}

zoom3DCore <- function(playState, foo)
{
    ## work out zoom factor by size of drag
    ## ideally this would set xlim/ylim/zlim to drag region?
    xlim <- rawXLim(playState)
    ylim <- rawYLim(playState)
    xfactor <- abs(diff(foo$coords$x)) / abs(diff(xlim))
    yfactor <- abs(diff(foo$coords$y)) / abs(diff(ylim))
    zoomfactor <- 1 / max(xfactor, yfactor)
    zoom <- callArg(playState, "zoom")
    if (is.null(zoom)) zoom <- 0.8
    callArg(playState, "zoom") <- signif(zoom * zoomfactor, 3)
    playReplot(playState)
}

zoomout3DCore <- function(playState, foo)
{
    zoom <- callArg(playState, "zoom")
    if (is.null(zoom)) zoom <- 0.8
    callArg(playState, "zoom") <- signif(zoom / 1.25, 3)
    playReplot(playState)
}

rotate3DCore <- function(playState, foo)
{
    ## work out current viewpoint (rotation)
    screen <- callArg(playState, "screen")
    if (is.null(screen)) screen <- list(z=40, x=-60)
    R.mat <- callArg(playState, "R.mat")
    if (is.null(R.mat)) R.mat <- diag(4)
    ## incorporate existing 'screen' into existing 'R.mat'
    R.mat <- ltransform3dMatrix(screen, R.mat)
    ## apply rotation defined by drag (direction and length)
    ## drag down corresponds to a positive 'x' arg
    ## drag right corresponds to a positive 'y' arg
    ## drag anticlockwise (in corner) corresponds to a positive 'z' arg
    pnl.x <- rawXLim(playState)
    pnl.y <- rawYLim(playState)
    x <- foo$coords$x
    y <- foo$coords$y
    xdelta <- diff(x)
    ydelta <- diff(y)
    angle <- atan2(y[2], x[2]) - atan2(y[1], x[1])
    dist <- c(max(abs(x/pnl.x)[1], abs(y/pnl.y)[1]),
              max(abs(x/pnl.x)[2], abs(y/pnl.y)[2]))
    ## TODO: avoid threshold for changing behaviour -- should be gradual
    if ((abs(angle) < pi/2) && all(dist > 0.5)) {
        rot <- list(z = 180 * angle / (2*pi))
    } else {
        ## TODO: should normalise by panel limits?
        rot <- list(y = xdelta * 180,
                    x = -ydelta * 180)
    }
    R.mat <- ltransform3dMatrix(rot, R.mat)
    R.mat <- round(R.mat, digits = 3)
    callArg(playState, "R.mat") <- call("matrix", c(R.mat), nc = 4)
    callArg(playState, "screen") <- list() ## replace default
    ## keep 3D scales on the same axes
    ## (it is confusing if they switch while rotating)
    ## ## no, bad for long axis labels
    #callArg(playState, "scpos") <- list(x = 1, y = 8, z = 4)
    playReplot(playState)
}

contextCore <- function(playState, foo, event)
{
    ## pop up context menu
  #cMenu <- gtkMenu()
  #cMenu$popup(button = event$button, activate.time = event$time)
  #cMenu["visible"] <- FALSE
    ## fill in menu items
    space <- foo$space
    canIdent <- playState$tmp$identify.ok
    if ((space != "page") && isTRUE(canIdent)) {
        foo$is.click <- TRUE
        foo <- playSelectData(playState, foo = foo, multiview = FALSE)
        id <- foo$subscripts
        pos <- foo$pos
        if (length(id) > 0) {
            id <- id[1]
            ## clicked on a data point, don't show general stuff
            cMenu <- gtkMenu()
            cMenu$popup(button = event$button, activate.time = event$time)
            #cMenu["visible"] <- TRUE
            ## action to add label to plot (current label value only)
            if (length(playState$labels) >= id) {
                item <- gtkMenuItem("Add label to plot:")
                item["sensitive"] <- FALSE
                cMenu$append(item)
                label <- toString(playState$labels[[id]])
                item <- gtkMenuItem(label)
                gSignalConnect(item, "activate",
                               function(widget, ...) {
                                   ## store newly identified points in playState
                                   playSetIDs(playState, id, pos = pos,
                                              type = "labelled",
                                              space = space,
                                              add = TRUE)
                               })
                cMenu$append(item)
            }
            ## covariate values
            cMenu$append(gtkSeparatorMenuItem())
            item <- gtkMenuItem("Annotate this data point:")
            item["sensitive"] <- FALSE
            cMenu$append(item)
            ## show values of x / y / other variables
            isLatt3D <-
                (playState$is.lattice &&
                 !is.null(playState$trellis$panel.args.common$scales.3d))
            x <- foo$x
            y <- foo$y
            if (is.numeric(x)) x <- round(x, 7)
            if (is.numeric(y)) y <- round(y, 7)
            dat <- getDataArg(playState)
            if (!is.null(dat)) {
                rn <- row.names(dat)
                itemtxt <- character()
                if (!is.null(rn)) {
                    itemtxt <- rn[id]
                }
                cn <- colnames(dat)
                if (length(cn) > 0) {
                    itemtxt <- c(itemtxt, paste(cn, ": ", dat[id,], sep = ""))
                }
                itemtxt <- sapply(itemtxt, toString, width = 33)
                for (txt in itemtxt) {
                    item <- gtkMenuItem(txt)
                    gSignalConnect(item, "activate",
                               function(widget, txt) {
                                   if (isLatt3D) {
                                       if (!gconfirm(paste("Annotation will be positioned in 2D only,",
                                                          "so will not match if you rotate the plot",
                                                          "(use 'Add label to plot', with 'Set labels to...' if you want that).",
                                                          "Continue?")))
                                           return()
                                   }
                                   ## add corresponding label
                                   annot <- call("panel.usertext", x, y, txt, pos = pos)
                                   playAnnotate(playState, annot, space = space)
                               }, data = txt)
                    cMenu$append(item)
                }
            } else {
                item <- gtkMenuItem(paste("x:", toString(x, width = 30)))
                cMenu$append(item)
                item <- gtkMenuItem(paste("y:", toString(y, width = 30)))
                cMenu$append(item)
            }
            ## "set labels to..."
            cMenu$append(gtkSeparatorMenuItem())
            aGroup <- playState$actionGroups[["PlotActions"]]
            cMenu$append(aGroup$getAction("SetLabelsTo")$createMenuItem())
            cMenu$append(gtkSeparatorMenuItem())
            cMenu$append(aGroup$getAction("SetLabelStyle")$createMenuItem())
            while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
            return()
        }
    }
    ## did not click on a point, so show general stuff
    #cMenu$destroy()
    cMenu <- playState$uiManager$getWidget("/ContextMenu")
    cMenu$popup(button = event$button, activate.time = event$time)
    while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
}
