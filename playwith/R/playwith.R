## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

## LICENSE
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version. See the file gpl-license.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


playwith <-
    function(expr,
             new = playwith.getOption("new"),
             title = NULL,
             labels = NULL,
             data.points = NULL,
             viewport = NULL,
             parameters = list(),
             tools = list(),
             init.actions = list(),
             preplot.actions = list(),
             update.actions = list(),
             ...,
             width = playwith.getOption("width"),
             height = playwith.getOption("height"),
             pointsize = playwith.getOption("pointsize"),
             eval.args = playwith.getOption("eval.args"),
             on.close = playwith.getOption("on.close"),
             modal = FALSE,
             link.to = NULL,
             playState = if (!new) playDevCur(),
             plot.call,
             main.function)
{
    if (!missing(plot.call) && !missing(expr))
        stop("give only one of 'expr' and 'plot.call'")
    if (missing(plot.call) && missing(expr)) {
        ## basic device mode
        expr <- quote({})
        if (missing(title))
            title <- "playwith (basic device mode)"
        show.call <- FALSE
    }
    if (missing(plot.call)) plot.call <- substitute(expr)
    if (is.expression(plot.call)) {
        plot.call <- if (length(plot.call) > 1)
            as.call(c(as.symbol("{"), plot.call)) else plot.call[[1]]
    }
    ## lattice calls can fall back to an update() wrapper
    if (is.symbol(plot.call) &&
        inherits(eval.parent(plot.call), "trellis"))
    {
        plot.call <- call("update", plot.call)
    }
    if (missing(main.function)) main.function <- NULL
    main.function <- substitute(main.function)
    if (is.language(main.function))
        main.function <- toString(main.function)
    ## check types
    if (!is.call(plot.call))
        stop("'expr' / 'plot.call' should be a call")
    if (!is.null(title) && !is.character(title)) stop("invalid 'title'")
    if (!is.null(viewport) && !is.list(viewport))
        viewport <- list(plot=viewport)
    if (!is.null(playState) && !inherits(playState, "playState"))
        stop("invalid 'playState'")
    if (!is.null(link.to) && !inherits(link.to, "playState"))
        stop("invalid 'link.to'")
    if (!is.list(init.actions))
        init.actions <- list(init.actions)
    if (!is.list(preplot.actions))
        preplot.actions <- list(preplot.actions)
    if (!is.list(update.actions))
        update.actions <- list(update.actions)
    ## playState is the <environment> encapsulating the plot, window and device
    cleanupStateEnv()
    if (!is.null(playState)) {
        stopifnot(inherits(playState, "playState"))
        if (is.null(title)) title <- playState$title
    }
    if (is.null(playState) || isTRUE(playState$keep)) {
        playState <- new.env(parent = emptyenv())
        class(playState) <- c("playState", "environment")
        ID <- basename(tempfile())
        playState$ID <- ID
        StateEnv[[ID]] <- playState
    }
    playState$tmp$plot.ready <- FALSE
    StateEnv$.current <- playState
    ## env is the <environment> containing local cached objects
    env <- new.env(parent = globalenv())
    ## work out evaluation rules
    envir <- parent.frame() ## where to look for variables
    invert.match <- FALSE
    if (is.list(eval.args)) {
        if (!is.null(eval.args$envir)) envir <- eval.args$envir
        if (!is.null(eval.args$invert)) invert.match <- eval.args$invert
        eval.args <- eval.args[[1]]
    }
    evalGlobals <- !is.na(eval.args)
    if (is.na(eval.args))
        eval.args <- (environmentName(envir) != "R_GlobalEnv")
    if (!identical(eval.args, FALSE)) {
        try(copyLocalArgs(plot.call, envir=envir, newEnv=env,
                          evalGlobals=evalGlobals, pattern=eval.args,
                          invert.match=invert.match))
    }
    ## check whether the window already exists
    myWin <- playState$win
    if (!is.null(myWin) && inherits(myWin, "GtkWindow")) {
        daAlloc <- playState$widgets$drawingArea$getAllocation()$allocation
        if (missing(width)) width <- daAlloc$width / 96
        if (missing(height)) height <- daAlloc$height / 96
        ## remove everything
        playState$tmp$devoff <- TRUE ## to avoid trigger close
        myWin$getChild()$destroy()
        myWin$present()
        myWin$resize(width * 96, height * 96)
    } else {
        ## create a new window
        myWin <- gtkWindow(show=FALSE)
        if (!inherits(myWin, "GtkWindow"))
            stop(paste("Could not create the GTK window.",
                       "Make sure you have recent versions of",
                       "RGtk2 and the GTK+ libraries."))
        ## set approx window size; NOTE: device size is adjusted below
        myWin["default-width"] <- width * 96
        myWin["default-height"] <- height * 96
        myWin["modal"] <- modal
        ## switch to GTK event loop while the window is in focus (for tooltips)
        myWin$addEvents(GdkEventMask["focus-change-mask"])
        gSignalConnect(myWin, "focus-in-event", gtkmain_handler,
                       data=playState)
        gSignalConnect(myWin, "focus-out-event", gtkmainquit_handler,
                       data=playState)
        gSignalConnect(myWin, "delete-event", gtkmainquit_handler,
                       data=playState)
        ## run user-defined close action
        gSignalConnect(myWin, "delete-event", window.close_handler,
                       data=playState)
    }
    if (!is.null(title)) myWin["title"] <- title
    myVBox <- gtkVBox()
    myWin$add(myVBox)
    playState$win <- myWin
    ## pass custom tools to UI manager
    optTools <- eval(playwith.getOption("custom.tools"))
    if (is.character(optTools)) optTools <- get(optTools)
    if (is.function(optTools)) optTools <- optTools(playState)
    cTools <- c(tools, optTools)
    ## extract any update / init actions attached to tools
    for (i in seq_along(cTools)) {
        xIni <- cTools[[i]]$init.action
        xPre <- cTools[[i]]$preplot.action
        xUpd <- cTools[[i]]$update.action
        if (!is.null(xIni)) {
            init.actions <- c(init.actions, xIni)
            cTools[[i]]$init.action <- NULL
        }
        if (!is.null(xPre)) {
            preplot.actions <- c(preplot.actions, xPre)
            cTools[[i]]$preplot.action <- NULL
        }
        if (!is.null(xUpd)) {
            update.actions <- c(update.actions, xUpd)
            cTools[[i]]$update.action <- NULL
        }
    }
    ## UI manager
    uiManager <- constructUIManager(playState, cTools)
    actionGroups <- uiManager$getActionGroups()
    names(actionGroups) <- sapply(actionGroups, gtkActionGroupGetName)
    ## construct menus
    menubar <- uiManager$getWidget("/MenuBar")
    menubar$show() ## location, layout behavior, padding
    menubar["visible"] <-
        isTRUE(playwith.getOption("show.menubar"))
    myVBox$packStart(menubar, expand=FALSE)
    ## construct the call toolbar
    callToolbar <- uiManager$getWidget("/CallToolbar")
    callToolbar["visible"] <-
        isTRUE(playwith.getOption("show.calltoolbar"))
    #callToolbar$setTooltips(TRUE) #playwith.getOption("show.tooltips"))
    callToolbar["toolbar-style"] <- GtkToolbarStyle["icons"]
    callToolbar["show-arrow"] <- FALSE
    ## merge in the address bar
    callEntry <- gtkComboBoxEntryNewText()
    ## load session history
    for (ihist in c(.PlaywithEnv$history, playState$history))
        callEntry$prependText(ihist)
    callEntry$show()
    ## "changed" emitted on typing and selection
    gSignalConnect(callEntry, "changed",
                   function(widget, playState)
                     if (widget["active"] > -1)
                       edit.call.inline_handler(widget$getChild(), playState),
                   data=playState)
    gSignalConnect(callEntry$getChild(), "activate",
                   edit.call.inline_handler, data=playState)
    item <- gtkToolItem()
    item$add(callEntry)
    item$setExpand(TRUE)
    callToolbar$insert(item, -1)
    callEditButton <- gtkButton(label="Edit call...")
    gSignalConnect(callEditButton, "clicked",
                   edit.call_handler, data=playState)
    item <- gtkToolItem()
    item$add(callEditButton)
    callToolbar$insert(item, -1)
    callToolbar$show()
    tbStyle <- GtkToolbarStyle[playwith.getOption("toolbar.style")]
    ## create the top toolbar
    topToolbar <- uiManager$getWidget("/TopToolbar")
#    topToolbar$setTooltips(TRUE)
    topToolbar["toolbar-style"] <- tbStyle
    ## create the bottom toolbar
    bottomToolbar <- uiManager$getWidget("/BottomToolbar")
#    bottomToolbar$setTooltips(TRUE)
    bottomToolbar["toolbar-style"] <- tbStyle
    ## create the left toolbar
    leftToolbar <- uiManager$getWidget("/LeftToolbar")
#    leftToolbar$setTooltips(TRUE)
    leftToolbar["toolbar-style"] <- tbStyle
    leftToolbar["orientation"] <- GtkOrientation["vertical"]
    ## create the right toolbar
    rightToolbar <- uiManager$getWidget("/RightToolbar")
#    rightToolbar$setTooltips(TRUE)
    rightToolbar["toolbar-style"] <- tbStyle
    rightToolbar["orientation"] <- GtkOrientation["vertical"]
    ## create the statusbar and coords readout
    statusbarBox <- gtkHBox()
    coordsLabel <- gtkLabel()
    coordsLabel["single-line-mode"] <- TRUE
    coordsLabel["selectable"] <- TRUE
    statusbarBox$packStart(coordsLabel, expand=FALSE)
    statusbar <- gtkStatusbar()
    statusbarBox$packStart(statusbar)
    statusbarBox["visible"] <-
        isTRUE(playwith.getOption("show.statusbar"))
    ## place toolbars in the window layout
    myVBox$packStart(callToolbar, expand=FALSE)
    myVBox$packStart(topToolbar, expand=FALSE)
    myHBox <- gtkHBox()
    myVBox$packStart(myHBox)
    myHBox$packStart(leftToolbar, expand=FALSE)
    myHBox$packEnd(rightToolbar, expand=FALSE)
    myVBox$packEnd(statusbarBox, expand=FALSE)
    myVBox$packEnd(bottomToolbar, expand=FALSE)
    ## create the plot area
    myDA <- gtkDrawingArea()
    myDA$addEvents(GdkEventMask["enter-notify-mask"]
                   + GdkEventMask["button-press-mask"]
                   + GdkEventMask["button-release-mask"]
                   + GdkEventMask["exposure-mask"])
    myHBox$packStart(myDA)
    myWin$show()
    asCairoDevice(myDA, pointsize = pointsize)
    ## note, this constraint is removed below
    dpi <- dev.size("px") / dev.size("in")
    myDA$setSizeRequest(width * dpi[1], height * dpi[2])
    ## try to force redraw
    gdkWindowProcessAllUpdates()
    while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
    ## need to regenerate coord spaces after resize
    gSignalConnect(myDA, "configure-event", configure_handler,
                   data=playState)
    gSignalConnect(myDA, "enter-notify-event", auto.reconfig_handler,
                   data=playState)
    gSignalConnect(myHBox, "remove", devoff_handler,
                   data=playState, after=TRUE)
    ## create the page scrollbar
    pageScrollBox <- gtkVBox(show=FALSE)
    myHBox$packStart(pageScrollBox, expand=FALSE)
    pageScrollBox$packStart(gtkLabel("Page"), expand=FALSE)
    pageEntry <- gtkEntry()
    pageEntry["width-chars"] <- 2
    gSignalConnect(pageEntry, "activate",
                   function(widget, playState) {
                       if (!isTRUE(playState$tmp$plot.ready)) return()
                       newPage <- round(as.numeric(widget["text"]))
                       if (newPage == playState$page) return()
                       playState$page <- newPage
                       playReplot(playState)
                   },
                   data=playState)
    pageScrollBox$packStart(pageEntry, expand=FALSE)
    pageScrollbar <- gtkVScrollbar()
    pageScrollbar["adjustment"] <-
        gtkAdjustment(value=1, lower=1, upper=1+1,
                      step.incr=1, page.incr=1, page.size=1)
    pageScrollbar["update-policy"] <- GtkUpdateType["delayed"]
    gSignalConnect(pageScrollbar, "value-changed",
                   function(widget, playState) {
                       if (!isTRUE(playState$tmp$plot.ready)) return()
                       newPage <- round(widget$getValue())
                       if (newPage == playState$page) return()
                       playState$page <- newPage
                       playReplot(playState)
                   },
                   data=playState)
    pageScrollBox$packStart(pageScrollbar)
    ## create the time/index scrollbar
    timeScrollBox <- gtkHBox(show=FALSE)
    myVBox$packStart(timeScrollBox, expand=FALSE)
    timeScrollBox$packStart(gtkLabel("Time"), expand=FALSE)
    timeEntry <- gtkEntry()
    timeEntry["width-chars"] <- 30
    gSignalConnect(timeEntry, "activate",
                   time.mode_entry_handler, data=playState)
    timeScrollBox$packStart(timeEntry, expand=FALSE)
    timeScrollbar <- gtkHScrollbar()
    timeScrollbar["adjustment"] <- gtkAdjustment()
    timeScrollbar["update-policy"] <- GtkUpdateType["delayed"]
    gSignalConnect(timeScrollbar, "value-changed",
                   time.mode_scrollbar_handler, data=playState)
    timeScrollBox$packStart(timeScrollbar)
    myHBox["resize-mode"] <- GtkResizeMode["queue"] ## does nothing?
    ## after resize, remove minimum size constraint from device
    myDA$setSizeRequest(-1, -1)
    ## store the state of this plot window in a new environment
    ## set per-window options -- can be replaced by explicit arguments
    playState$page <- 1
    playState$pages <- 1
    playState$is.lattice <- FALSE
    playState$is.ggplot <- FALSE
    playState$click.mode <- playwith.getOption("click.mode")
    playState$time.mode <- playwith.getOption("time.mode")
    playState$show.tooltips <- playwith.getOption("show.tooltips")
    playState$show.toolbars <- playwith.getOption("show.toolbars")
    playState$show.statusbar <- playwith.getOption("show.statusbar")
    playState$page.annotation <- playwith.getOption("page.annotation")
    playState$clip.annotations <- playwith.getOption("clip.annotations")
    playState$label.offset <- playwith.getOption("label.offset")
    playState$arrow <- playwith.getOption("arrow")
    ## store extra arguments (...) in the state object (playState)
    dots <- list(...)
    for (arg in names(dots)) {
        if (arg == "") next
        playState[[arg]] <- dots[[arg]]
    }
    ## initialise time.vector stuff
    if (!is.null(playState$time.vector)) {
        ## time.mode defaults to TRUE in this case
        if (is.null(dots$time.mode))
            playState$time.mode <- TRUE
        ## set current state variables
        env$cur.index <-
            if (!is.null(playState$cur.index)) {
                playState$cur.index
            } else if (!is.null(playState$cur.time)) {
                max(1, findInterval(playState$cur.time,
                                    playState$time.vector))
            } else 1
        env$cur.time <- playState$time.vector[env$cur.index]
        env$time.vector <- playState$time.vector
    }
    ## construct the state object (playState)
    playState$win <- myWin
    playState$dev <- dev.cur()
    playState$call <- plot.call
    playState$env <- env
    playState$tmp <- list()
    playState$labels <- labels
    playState$data.points <- data.points
    playState$viewport <- viewport
    playState$parameters <- parameters
    playState$tools <- tools
    playState$init.actions <- init.actions
    playState$preplot.actions <- preplot.actions
    playState$update.actions <- update.actions
    playState$on.close <- on.close
    playState$main.function <- main.function
    playState$pointsize <- pointsize
    playState$.args <-
        list(labels = labels,
             title = title)
    ## extras drawn on top of the plot
    playState$ids <- list()
    playState$annotations <- list()
    if (!is.null(link.to)) {
        playState$linked <- link.to$linked
        playState$linked$subscribers <-
            c(playState$linked$subscribers, playState)
        ## set click.mode from linked plot
        if (is.null(dots$click.mode))
            playState$click.mode <- link.to$tmp$click.mode
    } else {
        playState$linked <- new.env(parent = baseenv())
        playState$linked$ids <- list()
        playState$linked$subscribers <- list(playState)
    }
    playState$tmp$undoStack <- list()
    ## graphical user interface
    playState$uiManager <- uiManager
    playState$actionGroups <- actionGroups
    playState$widgets <-
        list(drawingArea = myDA,
             topToolbar = topToolbar,
             leftToolbar = leftToolbar,
             bottomToolbar = bottomToolbar,
             rightToolbar = rightToolbar,
             callToolbar = callToolbar,
             callEntry = callEntry,
             pageEntry = pageEntry,
             pageScrollbar = pageScrollbar,
             pageScrollBox = pageScrollBox,
             timeEntry = timeEntry,
             timeScrollbar = timeScrollbar,
             timeScrollBox = timeScrollBox,
             coordsLabel = coordsLabel,
             statusbar = statusbar,
             statusbarBox = statusbarBox,
             vbox = myVBox,
             hbox = myHBox)
    ## set initial values of any parameters
    for (i in seq_along(parameters)) {
        parname <- names(parameters)[i]
        parval <- parameters[[i]]
        if (is.list(parval)) parval <- parval[[1]]
        if (is.function(parval)) next
        assign(parname, parval[1], envir=env)
    }
    ## make dynamic parameter tools
    paramTbarNm <- paste("/", playwith.getOption("parameters.toolbar"), sep="")
    paramTbar <- playState$uiManager$getWidget(paramTbarNm)
    horiz <- (paramTbar["orientation"] == GtkOrientation["horizontal"])
    for (i in seq_along(parameters)) {
        parname <- names(parameters)[i]
        parval <- parameters[[i]]
        newTool <- try(parameterControlTool(playState, name=parname,
                                            value=parval, horizontal=horiz))
        if (inherits(newTool, "try-error")) next
        paramTbar$insert(newTool, -1)
    }
    ## try to force redraw
    gdkWindowProcessAllUpdates()
    while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
    ## fix toolbar sizes and hide empty toolbars
    blockRedraws({
        for (side in c("Top", "Left", "Bottom", "Right")) {
            nm <- paste("/", side, "Toolbar", sep = "")
            tbar <- playState$uiManager$getWidget(nm)
            if (isTRUE(playState$show.toolbars) &&
                length(tbar$getChildren()))
            {
                ## fix toolbar size
                sz <- tbar$getAllocation()$allocation
                if (side %in% c("Left", "Right")) {
                    tbar$setSizeRequest(sz$width, -1)
                } else {
                    tbar$setSizeRequest(-1, sz$height)
                }
            } else {
                tbar$hide()
            }
        }
    })
    ## try to force redraw
    gdkWindowProcessAllUpdates()
    while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
    ## do the plot
    playNewPlot(playState)
    invisible(playState)
}

playNewPlot <- function(playState = playDevCur())
    playReplot(playState, isNewPlot = TRUE)

playReplot <- function(playState = playDevCur(), isNewPlot = FALSE)
{
    if (isTRUE(playwith.getOption("catch.errors"))) {
        tryCatch(doPlayReplot(playState, isNewPlot = isNewPlot),
                 error = error_handler)
    } else {
        doPlayReplot(playState, isNewPlot = isNewPlot)
    }
}

doPlayReplot <- function(playState, isNewPlot = FALSE)
{
    if (!isNewPlot && isTRUE(playState$tmp$skip.redraws))
        return()
    playDevSet(playState)
    devAskNewPage(FALSE) ## 'ask' is always bad because of redraws
    playState$tmp$plot.ready <- FALSE
    on.exit(playState$tmp$plot.ready <- TRUE)
    grid.newpage()
    playPrompt(playState, NULL)
    ## disable toolbars until this is over
    playFreezeGUI(playState)
    on.exit(playThawGUI(playState))
    ## hide scrollbars if they are not needed
    if (!isTRUE(playState$time.mode)) {
        hideWidgetNoRedraw(playState, playState$widgets$timeScrollBox,
                           horiz = TRUE)
    }
    if (playState$pages == 1) {
        hideWidgetNoRedraw(playState, playState$widgets$pageScrollBox,
                           horiz = FALSE)
    }
    ## find which component of the call takes arguments (xlim etc)
    ## (also put main call into canonical form, with match.call)
    if (isNewPlot)
        updateMainCall(playState)
    ## pre-plot actions; note, can not assume $is.lattice etc is set
    preplotActions(playState)
    ## update address bar with current call
    updateAddressBar(playState)
    ## eval plot call
    ## (NOTE this will draw the plot UNLESS it is lattice or ggplot)
    result <- playState$result <-
        eval(playState$call, playState$env)
    ## detect lattice
    playState$is.lattice <- (inherits(result, "trellis"))
    playState$trellis <- if (playState$is.lattice) result
    ## detect ggplot
    playState$is.ggplot <- (inherits(result, "ggplot"))
    ## detect vcd
    playState$is.vcd <- (inherits(result, "structable"))
    ## detect base graphics
    ## TODO: use this elsewhere
    playState$is.base <- (!playState$is.lattice &&
                          !playState$is.ggplot &&
                          !playState$is.vcd &&
                          is.null(playState$viewport))
    ## lattice needs subscripts to correctly identify points.
    if (isNewPlot && playState$is.lattice &&
        prod(dim(playState$trellis)) > 1)
    {
        if (is.null(playState$trellis$panel.args[[1]]$subscripts)) {
            warning(paste("Call may need subscripts = TRUE",
                          "to correctly identify data points."))
        }
    }
    ## initialisation actions for a new plot (see uiManager.R)
    if (isNewPlot)
        initActions(playState)
    ## set back to this device, since may have switched during plot
    playDevSet(playState)
    playState$pages <- 1
    if (playState$is.lattice) {
        ## work out number of pages
        playState$pages <- npages(result)
        if (playState$page > playState$pages)
            playState$page <- 1
        ## plot trellis object (specified page only)
        plotOnePage(result, page = playState$page)
        ## need to store this, it refers to last plot only!
        playState$tmp$currentLayout <-
            trellis.currentLayout(which="packet")
        ## use lattice style settings attached to trellis object
        ## for any annotations (in updateActions(), below)
        if (!is.null(result$par.settings)) {
            opar <- trellis.par.get()
            trellis.par.set(result$par.settings)
            on.exit(trellis.par.set(opar, strict = TRUE), add = TRUE)
        }
    }
    if (playState$is.ggplot) {
        ## plot ggplot object
        print(result)
        ## typically want: playState$viewport <- list(plot = "panel_1_1")
        vpNames <- grid.ls(viewports=TRUE, grobs=FALSE, print=FALSE)$name
        panelNames <- vpNames[grep("panel-", vpNames)]
        panelNames <- unique(panelNames)
        tmp.vp <- as.list(panelNames)
        names(tmp.vp) <- panelNames
        names(tmp.vp)[1] <- "plot"
        playState$viewport <- tmp.vp
    }
#    if (inherits(result, "structable")) {
#        ## can't really interact in this case
#        playState$viewport <- NULL
#    }
    ## store coordinate system(s)
    generateSpaces(playState)
    playState$tmp$plot.ready <- TRUE
    ## update toolitem and menuitem states (see uiManager.R)
    updateActions(playState)
    ## and update the pages scrollbar
    pages_post.plot.action(playState$widgets$pageScrollBox,
                           playState=playState)
    invisible(result)
}

updateAddressBar <- function(playState)
{
    ## add current call to "address bar" and set up associated tools
    callTxt <- ""
    if (object.size(playState$call) < 50000) {
        widg <- playState$widgets
        callTxt <- deparseOneLine(playState$call, control=
                                  playwith.getOption("deparse.options"))
        if (is.null(playState$.args$title)) playState$win["title"] <-
            toString(callTxt, width=34)
        oldCallTxt <- widg$callEntry$getActiveText()
        if ((widg$callEntry["active"] == -1) || (callTxt != oldCallTxt)) {
            ## a new call: edited inline OR playState$call modified
            widg$callEntry$prependText(callTxt)
            widg$callEntry["active"] <- 0
            ## record in history
            playState$history <- c(playState$history, callTxt)
            ## remove any later history (branching from a previous state)
            histLev <- playState$tmp$call.history.level
            if (any(histLev > 0)) {
                for (i in seq(histLev-1, 0)+1)
                  widg$callEntry$removeText(i)
                #playState$history <- head(playState$history, -histLev)
            }
        }
        playState$tmp$call.history.level <- widg$callEntry["active"]
    }
}

generateSpaces <- function(playState)
{
    ## enumerate spaces in the current plot
    ## (named list of viewports)
    playState$spaces <- list()
    upViewport(0)
    curVps <- grid.ls(grobs = FALSE, viewports = TRUE, print = FALSE)$name
    if (!is.null(playState$viewport)) {
        ## grid graphics plot
        playState$spaces <- names(playState$viewport)
    } else if (playState$is.lattice) {
        ## lattice plot
        packets <- playState$tmp$currentLayout
        playState$spaces <- paste("packet", packets[packets > 0])
    } else if (playState$is.base) {
        ## base graphics plot
        playState$spaces <- "plot"
        ## use gridBase to make viewports
        if (length(playState$tmp$baseVps$plot.clip.off)) {
            if ("plot.clip.off" %in% curVps) {
                downViewport("plot.clip.off")
                popViewport(0)
            }
        }
        ## suppress warnings about log scale
        vps <- suppressWarnings(baseViewports())
        vps$plot$name <- "plot"
        vps$plot$clip <- TRUE
        vps$plot.clip.off <-
            viewport(xscale=par("usr")[1:2],
                     yscale=par("usr")[3:4],
                     clip="off", name = "plot.clip.off")
        playState$tmp$baseVps <- vps
        pushViewport(do.call("vpStack", vps))
    }
    upViewport(0)
    ## create a top-level viewport with normalised coordinates
    ## yscale origin is at top, to be consistent with device coordinates
    if (("pageAnnotationVp" %in% curVps) == FALSE) {
        if (playState$is.lattice)
            downViewport(trellis.vpname("toplevel"))
        pushViewport(viewport(name = "pageAnnotationVp",
                              yscale = c(1, 0)))
        upViewport(0)
    }
    ## store coordinate transformations for each space
    playState$tmp$spaceLimDevice <- list()
    for (space in playState$spaces) {
        ## bounds in device coordinates (pixels)
        playState$tmp$spaceLimDevice[[space]] <-
            playDo(playState,
                   convertToDevicePixels(x = unit(0:1, "npc"),
                                         y = unit(0:1, "npc")),
                   space = space)
    }
    playState$tmp$need.reconfig <- FALSE
}

### Error handler

error_handler <- function(e)
{
    if (inherits(e, "error")) {
        callText <- toString(deparseOneLine(conditionCall(e)),
                             width = 200)
        msg <- paste("Error: ", conditionMessage(e),
                     "\n\nThe error occurred in: \n",
                     callText, sep="")
        gmessage.error(msg)
    }
    e
}

### Window signal handlers

pages_post.plot.action <- function(widget, playState)
{
    widg <- playState$widgets
    if (playState$pages > 1) {
        widg$pageScrollbar["adjustment"]["upper"] <- playState$pages+1
        widg$pageScrollbar["adjustment"]["value"] <- playState$page
        widg$pageEntry["text"] <- toString(playState$page)
        widg$pageScrollBox["sensitive"] <- TRUE
        widg$pageScrollbar["sensitive"] <- TRUE
        widg$pageEntry["sensitive"] <- TRUE
        if (widg$pageScrollBox["visible"] == FALSE) {
            blockRedraws(widg$pageScrollBox$show())
        }
    } else {
                                        #widg$pageScrollBox$hide()
        widg$pageScrollbar["sensitive"] <- FALSE
        widg$pageEntry["sensitive"] <- FALSE
    }
}

window.close_handler <- function(widget, event, playState)
{
    if (isTRUE(playState$keep)) {
        ans <- gconfirm("This plot is marked to keep open. Really close?",
                        parent = playState$win)
        if (!isTRUE(ans)) return(TRUE) ## do not close
    }
    if (!is.null(playState$on.close)) {
        foo <- try(playState$on.close(playState))
        ## if on.close() returns TRUE, do not close the window
        if (isTRUE(foo)) return(TRUE)
    }
## More annoying than useful I think:
#    if (length(playState$linked$subscribers) > 1) {
#        ans <- gconfirm("Also close linked plots?",
#                        parent = playState$win)
#        if (isTRUE(ans)) {
#            lapply(playState$linked$subscribers,
#                   playDevOff)
#            return(FALSE)
#        }
#    }
    ## close the window and clean up
    playDevOff(playState)
    return(FALSE)
}

configure_handler <- function(widget, event, playState)
{
    playState$tmp$need.reconfig <- TRUE
    return(FALSE)
}

auto.reconfig_handler <- function(widget, event, playState)
{
    ## avoid weird stack smash
    if (length(playState$is.lattice) == 0) return(FALSE)
    if (!isTRUE(playState$tmp$plot.ready)) return(FALSE)
    if (isBasicDeviceMode(playState)) {
        ## do not know when the plot is updated
        ## so need to keep regenerating data space
        playState$tmp$need.reconfig <- TRUE
    }
    if (playState$tmp$need.reconfig) {
        generateSpaces(playState)
    }
    return(FALSE)
}

devoff_handler <- function(widget, event, playState)
{
    ## this handles dev.off()
    ## destroy the window, but store a flag to avoid destroying twice
    if (isTRUE(playState$tmp$devoff)) return(FALSE)
    playState$tmp$devoff <- TRUE
    playDevOff(playState)
    return(FALSE)
}

gtkmain_handler <- function(widget, event, playState)
{
  if (!isTRUE(playState$show.tooltips))
    return(gtkmainquit_handler(widget, event, playState))
  ## switch to GTK event loop while the window is in focus (for tooltips)
  if (!isTRUE(playState$tmp$gtkMain)) {
    playState$tmp$gtkMain <- TRUE
    gtkMain()
  }
  return(FALSE)
}

gtkmainquit_handler <- function(widget, event, playState)
{
  if (isTRUE(playState$tmp$gtkMain)) {
    playState$tmp$gtkMain <- FALSE
    gtkMainQuit()
  }
  return(FALSE)
}

## General utility functions

## i'm pretty sure this behaviour used to be in base R:
toString.function <- function(x, ...)
    toString(deparse(x), ...)

deparseOneLine <-
    function(expr, width.cutoff = 500, ...)
{
    tmp <- deparse(expr, width.cutoff = width.cutoff, ...)
    indents <- attr(regexpr("^ *", tmp), "match.length")
    breaks <- c(diff(indents) <= 0, FALSE)
    tmp <- gsub("^ +", "", tmp)
    tmp <- gsub(" +$", "", tmp)
    breaks[c(tmp[-1]=="{", FALSE)] <- F
    tmp <- paste(tmp, ifelse(breaks, ";", ""), sep="", collapse=" ")
    tmp <- gsub("\\{;", "\\{", tmp)
    tmp <- gsub(";\\}", " \\}", tmp)
    tmp <- gsub(";\\{", " \\{", tmp)
    ## update: need this for long inline vectors:
    tmp <- gsub(";,", ",", tmp)
    tmp <- gsub(",;", ",", tmp)
    tmp
}

copyLocalArgs <-
    function(the.call,
             envir = parent.frame(),
             newEnv,
             evalGlobals = FALSE,
             pattern = TRUE,
             invert.match = FALSE)
{
    stopifnot(is.call(the.call) || is.list(the.call) || is.expression(the.call))
    isMatch <- !invert.match
    for (i in seq_along(the.call)) {
        this.arg <- the.call[[i]]
        if (is.call(the.call) && (i == 1)) {
            if (mode(this.arg) == "name") {
                callname <- as.character(this.arg)
                ## skip base operators, unlikely to be locally redefined
                if (exists(callname, baseenv(), mode="function"))
                    next
            } else {
                ## anonymous function call
                next
            }
        }

        if (is.call(this.arg)) {
            if (mode(this.arg[[1]]) == "name") {
                argcallname <- as.character(this.arg[[1]])
                ## skip function definitions
                ## (could skip function entirely, but it may refer
                ##  to variables in outer context)
                if (argcallname == "function")
                    this.arg <- this.arg[[3]]
                ## skip assigned variables
                if (argcallname == "<-")
                    this.arg <- this.arg[[3]]
                ## skip literal symbol in "$" extractor
                if (argcallname == "$")
                    this.arg <- this.arg[[2]]
            }
        }

        if (mode(this.arg) %in% c("call", "(", "list", "expression")) {
            ## call recursively...
            copyLocalArgs(this.arg, envir=envir, newEnv=newEnv,
                          evalGlobals=evalGlobals, pattern=pattern,
                          invert.match=invert.match)
        } else if (mode(this.arg) == "name") {
            this.name <- as.character(this.arg)
            ## check if this name already exists in local env
            if (exists(this.name, newEnv, inherits=F))
                next
            ## check that the name matches the pattern
            if (!isTRUE(pattern) &&
                (any(grep(pattern, this.name))) != isMatch)
                next
            ## check if name exists in 'envir' or its parents
            ## (up to global env only, i.e. local objects)
            testenv <- envir
            hit <- FALSE
            while (TRUE) {
                if (exists(this.name, testenv, inherits=FALSE)) {
                    assign(this.name, get(this.name, testenv), newEnv)
                    hit <- TRUE
                    break
                }
                testenv <- parent.env(testenv)
                if (environmentName(testenv) %in%
                    c("R_GlobalEnv", "R_EmptyEnv", "base"))
                    break
            }
            if (hit) next
            ## eval globals
            if (evalGlobals &&
                exists(this.name, globalenv(), inherits=FALSE))
                assign(this.name, get(this.name, testenv), newEnv)
        }
        ## leave constants alone
    }
}

npages <- function(x) {
    stopifnot(inherits(x, "trellis"))
    ## work out number of pages that would be plotted
    ## to display trellis object 'x'
    nPackets <- prod(dim(x))
    ## by default, first two dimensions
    ## (conditioning variables) shown on each page
    nPanels <- prod(head(dim(x), 2))
    ## but if an explicit 'layout' is given...
    if (!is.null(x$layout)) {
        nPanels <- x$layout[1] * x$layout[2]
        if (x$layout[1] == 0) nPanels <- x$layout[2]
    }
    ## TODO: what about 'skip'?
    nPages <- ceiling(nPackets / nPanels)
    nPages
}

plotOnePage <- function(x, page, ...)
{
    stopifnot(inherits(x, "trellis"))
    n <- page
    if (is.null(x$layout)) {
        if (length(dim(x)) > 2)
            x$layout <- dim(x)[1:2]
    }
    if (!is.null(x$layout))
        x$layout[3] <- 1
    ## based on code by Deepayan Sarkar
    packet.panel.pageN <- function(..., page)
        packet.panel.default(..., page = page + n - 1)
    plot(x, packet.panel = packet.panel.pageN, ...)
}

## unused
recursive.as.list.call <- function(x) {
    stopifnot(is.call(x))
    x <- as.list(x)
    lapply(x, function(z) if (is.call(z))
           recursive.as.list.call(z) else z)
}
