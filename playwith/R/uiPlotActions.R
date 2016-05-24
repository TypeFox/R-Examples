## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

plotActionGroup <- function(playState)
{
    entries <-
        list( ## : name, stock icon, label, accelerator, tooltip, callback
             list("PlotSettings", "gtk-preferences", "Plot _settings", "<Ctrl>semicolon", "Change the plot type and settings", plot.settings_handler),
             list("Zoomin", "gtk-zoom-in", "Zoom _in", "<Ctrl>bracketright", "Show half the scale of x and y", zoomin_handler),
             list("Zoomout", "gtk-zoom-out", "Zoom _out", "<Ctrl>bracketleft", "Show twice the scale of x and y", zoomout_handler),
             list("Zoomfit", "gtk-zoom-fit", "_Reset scales", "<Ctrl>0", "Revert to default plot region", zoomfit_handler),
             list("ZeroY", "gtk-goto-bottom", "Full y scale", "<Ctrl>Return", "Show the full y (response) scale starting from zero", zero.y_handler),
             list("ZeroX", "gtk-goto-first", "Full x scale", "<Ctrl>BackSpace", "Show the full x (domain) scale starting from zero", zero.x_handler),
             ## identify (uiIdentifyActions.R)
             list("SetLabelsTo", "gtk-index", "Set _labels to...", "<Ctrl>L", NULL, set.labels_handler),
             list("IdTable", "gtk-index", "Select from _table...", "<Ctrl>J", "Select points from a table", id.table_handler),
             list("FindLabels", "gtk-find", "_Find...", "<Ctrl>F", "Find points with labels matching...", id.find_handler),
             list("SaveIDs", NULL, "_Save IDs...", NULL, "_Save current IDs to an object", save.ids_handler),
             ## annotations (uiAnnotationActions.R)
             list("Legend", "gtk-sort-ascending", "Legend", NULL, "Place a legend", legend_handler),
             list("EditAnnotations", "gtk-edit", "Edit annot.", "<Ctrl><Shift>E", "Edit annotations (including arrows) code", edit.annotations_handler),
             list("UndoAnnotation", "gtk-undo", "Undo annot.", "<Ctrl>Z", "Remove last annotation", undo.annotation_handler),
             list("Clear", "gtk-clear", "Clear", "<Shift>Delete", "Remove labels and annotations", clear_handler),
             ## style (uiStyleActions.R)
             list("SetPointLineStyle", NULL, "Set _point/line style...", "<Alt>1", NULL, set.point.line.style_handler),
             list("SetBrushStyle", NULL, "Set _brush style...", "<Alt>2", NULL, set.brush.style_handler),
             list("SetLabelStyle", NULL, "Set _label style...", "<Alt>3", NULL, set.label.style_handler),
             list("SetArrowStyle", NULL, "Set _arrow style...", "<Alt>4", NULL, set.arrow.style_handler),
             list("DefaultTheme", NULL, "Default", NULL, NULL, set.default.theme_handler),
             ## grobs (uiGrobActions.R)
             list("GrobInspector", "gtk-properties", "_Grob inspector", "<Ctrl>question", NULL, grob.inspector_handler)
             )

    toggleEntries <-
        list( ## : name, stock icon, label, accelerator, tooltip, callback, active
             list("Expand", "gtk-fullscreen", "_Panel", NULL, "Choose a panel to expand to fill the figure (for further interaction)", expand_handler, FALSE),
             ## style shortcuts (uiStyleActions.R)
             list("StyleSolidPoints", NULL, "Solid points", NULL, NULL, style.solid.points_handler, FALSE),
             list("StyleTransPoints", NULL, "Translucent points", NULL, NULL, style.trans.points_handler, FALSE),
             list("StyleThickLines", NULL, "Thick lines", NULL, NULL, style.thick.lines_handler, FALSE),
             ## options (uiOptionsActions.R)
             list("TimeMode", "gtk-media-forward-ltr", "_Time mode", "<Ctrl>T", "Time mode: scroll along the x axis", time.mode_handler, FALSE)
             )

    ## see uiClickActions.R
    clickModeEntries <-
        list( ## : name, stock icon, label, accelerator, tooltip, value
             list("Zoom", "gtk-zoom-in", "_Navigate", "<Ctrl>space", "Zoom in and out or rotate (3D)", 0),
             list("Pan", "gtk-jump-to", "_Pan", "<Alt>space", "Pan (move the plot area)", 1),
             list("Identify", "gtk-info", "_Identify", "<Ctrl>I", "Identify data points (add labels)", 2),
             list("Brush", "gtk-media-record", "_Brush", "<Ctrl>B", "Brush (highlight) data points", 3),
             list("Annotation", "gtk-italic", "_Annotate", "<Ctrl>apostrophe", "Add text to the plot", 4),
             list("Arrow", "gtk-connect", "Arro_w", "<Ctrl>bar", "Add arrows to the plot", 5),
             list("Line", NULL, "_Line", "<Alt>bar", "Add lines to the plot", 6),
             list("Rect", NULL, "_Rect", "<Alt>apostrophe", "Add rectangles to the plot", 7)
             )

    ## construct action group with playState passed to callbacks
    aGroup <- gtkActionGroupNew("PlotActions")
    aGroup$addActions(entries, playState)
    aGroup$addToggleActions(toggleEntries, playState)
    aGroup$addRadioActions(clickModeEntries, 0,
                           on.change = clickmode.change_handler,
                           playState)
    aGroup
}

clickModeValues <- function() {
    list("Zoom" = 0,
         "Pan" = 1,
         "Identify" = 2,
         "Brush" = 3,
         "Annotation" = 4,
         "Arrow" = 5,
         "Line" = 6,
         "Rect" = 7)
}

updatePlotActions <- function(playState)
{
    aGroup <- playState$actionGroups[["PlotActions"]]
    hasArgs <- playState$accepts.arguments
    isLatt <- playState$is.lattice
    isLatt3D <-
        (isLatt &&
         !is.null(playState$trellis$panel.args.common$scales.3d))
    isVCD <- playState$is.vcd
    canNav <- (hasArgs && !isVCD &&
               !(playState$callName %in%
                 c("splom", "marginal.plot")))
    ## PlotSettings
    aGroup$getAction("PlotSettings")$setSensitive(hasArgs)
    ## Zoom / Pan
    aGroup$getAction("Zoom")$setSensitive(canNav)
    aGroup$getAction("Pan")$setSensitive(canNav)
    ## Zoomfit
    nonFit <- FALSE
    if (canNav) {
        nonFit <-
            hasArgs && (!is.null(callArg(playState, "xlim")) ||
                        !is.null(callArg(playState, "ylim")))
        if (isLatt3D)
            nonFit <- (nonFit ||
                       !is.null(callArg(playState, "zlim")) ||
                       !is.null(callArg(playState, "zoom")) ||
                       !is.null(callArg(playState, "screen")) ||
                       !is.null(callArg(playState, "R.mat")))
    }
    aGroup$getAction("Zoomfit")$setSensitive(nonFit)
    aGroup$getAction("Zoomfit")$setVisible(nonFit)
    ## ZeroY
    eps <- .Machine$double.eps * 2
    nonZeroY <- FALSE
    if (canNav) {
        ylim <- rawYLim(playState)
        if (isLatt) {
            ylim <- playState$trellis$y.limits
            if (is.list(ylim)) ylim <- ylim[[1]]
            if (is.character(ylim)) ylim <- c(0, 0)
        }
        if (isLatt3D)
            ylim <- playState$trellis$panel.args.common$zlim
        nonZeroY <- (min(ylim) > eps) || (max(ylim) < -eps)
    }
    aGroup$getAction("ZeroY")$setSensitive(nonZeroY)
    aGroup$getAction("ZeroY")$setVisible(nonZeroY)
    ## ZeroX
    nonZeroX <- FALSE
    if (canNav) {
        xlim <- rawXLim(playState)
        if (isLatt) {
            xlim <- playState$trellis$x.limits
            if (is.list(xlim)) xlim <- xlim[[1]]
            if (is.character(xlim)) xlim <- c(0, 0)
        }
        if (isLatt3D)
            xlim <- playState$trellis$panel.args.common$xlim
        nonZeroX <- (min(xlim) > eps) || (max(xlim) < -eps)
    }
    aGroup$getAction("ZeroX")$setSensitive(nonZeroX)
    aGroup$getAction("ZeroX")$setVisible(nonZeroX)
    ## Expand
    hasPanels <- isLatt && (length(playState$tmp$currentLayout) > 1)
    expandActive <- aGroup$getAction("Expand")$getActive()
    aGroup$getAction("Expand")$setVisible(expandActive || hasPanels)
}

zoomout_handler <- function(widget, playState)
{
    isLatt <- playState$is.lattice
    isLatt3D <- isLatt && !is.null(playState$trellis$panel.args.common$scales.3d)
    if (isLatt3D) {
        zoom <- callArg(playState, "zoom")
        if (is.null(zoom)) zoom <- 0.8
        callArg(playState, "zoom") <- signif(zoom / 1.25, 3)
        playReplot(playState)
        return()
    }
    nav.x <- TRUE
    ## in time.mode, only zoom along x-axis
    nav.y <- (!isTRUE(playState$time.mode) &&
              is.null(playState$time.vector))
    ## find existing scales
    xlim <- rawXLim(playState)
    ylim <- rawYLim(playState)
    ## zoom out: make range twice the size
    if (nav.x) xlim <- xlim + diff(xlim) * c(-0.5, 0.5)
    if (nav.y) ylim <- ylim + diff(ylim) * c(-0.5, 0.5)
    ## this converts from raw numeric to original format (including unlog)
    if (nav.x) rawXLim(playState) <- xlim
    if (nav.y) rawYLim(playState) <- ylim
    playReplot(playState)
}

zoomin_handler <- function(widget, playState)
{
    isLatt <- playState$is.lattice
    isLatt3D <- isLatt && !is.null(playState$trellis$panel.args.common$scales.3d)
    if (isLatt3D) {
        zoom <- callArg(playState, "zoom")
        if (is.null(zoom)) zoom <- 0.8
        callArg(playState, "zoom") <- signif(zoom * 1.25, 3)
        playReplot(playState)
        return()
    }
    nav.x <- TRUE
    ## in time.mode, only zoom along x-axis
    nav.y <- (!isTRUE(playState$time.mode) &&
              is.null(playState$time.vector))
    ## find existing scales
    xlim <- rawXLim(playState)
    ylim <- rawYLim(playState)
    ## zoom in: make range half the size
    if (nav.x) xlim <- xlim + diff(xlim) * c(0.25, -0.25)
    if (nav.y) ylim <- ylim + diff(ylim) * c(0.25, -0.25)
    ## this converts from raw numeric to original format (including unlog)
    if (nav.x) rawXLim(playState) <- xlim
    if (nav.y) rawYLim(playState) <- ylim
    playReplot(playState)
}

zoomfit_handler <- function(widget, playState)
{
    isLatt <- playState$is.lattice
    isLatt3D <- isLatt && !is.null(playState$trellis$panel.args.common$scales.3d)
    ## update scales
    callArg(playState, "xlim") <- NULL
    callArg(playState, "ylim") <- NULL
    if (isLatt3D) {
        callArg(playState, "zlim") <- NULL
        callArg(playState, "zoom") <- NULL
        callArg(playState, "screen") <- NULL
        callArg(playState, "R.mat") <- NULL
    }
    playReplot(playState)
}

zero.x_handler <- function(widget, playState)
{
    if (playState$is.lattice) {
        is3D <- !is.null(playState$trellis$panel.args.common$scales.3d)
        if (is3D) {
            ## in 3D case, zero both domain variables ("x" and "y")
            xlim <- playState$trellis$panel.args.common$xlim
            ylim <- playState$trellis$panel.args.common$ylim
            xlim[which.min(abs(xlim))] <- 0
            ylim[which.min(abs(ylim))] <- 0
            callArg(playState, "xlim") <- signif(xlim, 7)
            callArg(playState, "ylim") <- signif(ylim, 7)
            playReplot(playState)
            return()
        } else {
            xlim <- playState$trellis$x.limits
            if (is.character(xlim)) return()
            if (is.list(xlim)) {
                if (!all(lapply(xlim, is.numeric)))
                    return()
                makeScalesArgAList(playState)
                xlim <- lapply(xlim, function(lim)
                           {
                               lim[which.min(abs(lim))] <- 0
                               signif(lim, 7)
                           })
                callArg(playState, quote(scales$x$limits)) <-
                    xlim
                playReplot(playState)
                return()
            }
        }
    } else {
        xlim <- rawXLim(playState)
    }
    xlim[which.min(abs(xlim))] <- 0
    rawXLim(playState) <- signif(xlim, 7)
    playReplot(playState)
}

zero.y_handler <- function(widget, playState)
{
    if (playState$is.lattice) {
        is3D <- !is.null(playState$trellis$panel.args.common$scales.3d)
        if (is3D) {
            ## in 3D case, zero the response variable ("z")
            zlim <- playState$trellis$panel.args.common$zlim
            zlim[which.min(abs(zlim))] <- 0
            callArg(playState, "zlim") <- signif(zlim, 7)
            playReplot(playState)
            return()
        } else {
            ylim <- playState$trellis$y.limits
            if (is.character(ylim)) return()
            if (is.list(ylim)) {
                if (!all(lapply(ylim, is.numeric)))
                    return()
                makeScalesArgAList(playState)
                ylim <- lapply(ylim, function(lim)
                           {
                               lim[which.min(abs(lim))] <- 0
                               signif(lim, 7)
                           })
                callArg(playState, quote(scales$y$limits)) <-
                    ylim
                playReplot(playState)
                return()
            }
        }
    } else {
        ylim <- rawYLim(playState)
    }
    ylim[which.min(abs(ylim))] <- 0
    rawYLim(playState) <- signif(ylim, 7)
    playReplot(playState)
}

expand_handler <- function(widget, playState)
{
    playDevSet(playState)
    ## check new expanded setting
    if (widget["active"]) {
        playPrompt(playState,
                   paste("Click on a panel to expand;",
                         "Right-click or Esc to cancel."))
        on.exit(playPrompt(playState, NULL))
        newFocus <- trellis.focus()
        if (is.null(newFocus) || all(newFocus == 0)) {
            widget["active"] <- FALSE
            return()
        }
        playState$tmp$old.call.layout <- callArg(playState, "layout")
        callArg(playState, "layout") <- c(0,1,1)
        playState$tmp$old.page <- playState$page
        playState$page <- packet.number()
        playState$tmp$old.pages <- playState$pages
    } else {
        if (is.null(playState$tmp$old.page)) return()
        callArg(playState, "layout") <- playState$tmp$old.call.layout
        playState$page <- playState$tmp$old.page
        playState$pages <- playState$tmp$old.pages
        playState$tmp$old.page <- NULL
    }
    playReplot(playState)
}
