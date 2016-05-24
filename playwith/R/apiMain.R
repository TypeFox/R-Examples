# playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer


playDevCur <- function()
    StateEnv$.current ## may be NULL

playDevList <- function()
{
    foo <- as.list(StateEnv)
    names(foo) <- lapply(foo, function(x)
                         toString(x$win["title"]))
    foo
}

playDevSet <- function(playState = playDevCur())
{
    stopifnot(inherits(playState, "playState"))
    StateEnv$.current <- playState
    playState$tmp$old.dev <- dev.cur()
    dev.set(playState$dev)
}

playDevOff <- function(playState = playDevCur())
{
    ## save local history to session history
    .PlaywithEnv$history <-
        c(.PlaywithEnv$history, playState$history)
    playState$history <- NULL
    ## TODO: should this run the close action?
    if (inherits(playState$win, "GtkWindow"))
        try(playState$win$destroy())#, silent=TRUE)
    ## it seems that memory is not freed! (R2.7.1)
    rm(list=ls(playState), envir=playState)
    cleanupStateEnv()
}

print.playState <- function(x, ...)
{
    stopifnot(inherits(x, "playState"))
    title <- "(invalid)"
    if (inherits(x$win, "GtkWindow"))
        title <- toString(x$win["title"])
    cat(paste("<playState: ", title, ">\n", sep=""))
}

cleanupStateEnv <- function()
{
    for (ID in ls(StateEnv)) {
        if (!inherits(StateEnv[[ID]], "playState")) next
        if (!inherits(StateEnv[[ID]]$win, "GtkWindow")) {
            ## window is defunct
            rm(list=ID, envir=StateEnv)
        }
    }
    ## select a new 'current' if it is invalid
    if (!inherits(StateEnv$.current$win, "GtkWindow")) {
        StateEnv$.current <- if (length(ls(StateEnv)))
            StateEnv[[ ls(StateEnv)[1] ]] else NULL
    }
}

playwith.history <- function(max.show = 100, ...)
{
    txt <- 
        c(.PlaywithEnv$history,
          unlist(lapply(playDevList(), function(x) x$history)))
    if (length(txt) == 0) {
        message("No history to display.")
        return(invisible())
    }
    file2 <- tempfile("Rplaywithhist")
    inds <- tail(seq_along(txt), max.show)
    writeLines(txt[inds], file2)
    file.show(file2, title = "playwith history", delete.file = TRUE)
}

callArg <- function(playState, arg, eval = TRUE, data = NULL)
{
    if (is.symbol(arg)) arg <- as.character(arg)
    getx <- if (is.numeric(arg)) paste('[[', arg+1, ']]', sep="")
    else if (is.character(arg)) paste('[["', arg, '", exact=TRUE]]', sep="")
    else paste("$", deparseOneLine(arg), sep="")
    mainCall <- mainCall(playState)
    zap <- eval(parse(text=paste("mainCall", getx, sep="")))
    if (eval == FALSE) return(zap)
    if (mode(zap) == "expression") return(zap)
    if (is.null(data))
        eval(zap, envir=playState$env, enclos=parent.frame())
    else
        eval(zap, envir=data, enclos=playState$env)
}

"callArg<-" <- function(playState, arg, value)
{
    if (is.symbol(arg)) arg <- as.character(arg)
    if (is.null(arg)) return()
    getx <- if (is.numeric(arg)) paste('[[', arg+1, ']]', sep="")
    else if (is.character(arg)) paste("$", arg, sep="")
    else paste("$", deparseOneLine(arg), sep="")
    mainCall <- mainCall(playState)
    zap <- parse(text=paste("mainCall", getx, sep=""))[[1]]
    zap <- call("<-", zap, quote(value))
    eval(zap, enclos=parent.frame())
    ## instantiate implicit lists as language objects
    ## this is required for e.g. lattice's scales$x$at <- quote(qnorm(...))
    ## easiest way is just to deparse without showAttributes and then parse
    if (is.language(arg)) {
        tmp <- try( parse(text = deparse(mainCall,
                           control = playwith.getOption("deparse.options")
                           ))[[1]] )
        if (!inherits(tmp, "try-error"))
            mainCall <- tmp
    }
    mainCall(playState) <- mainCall
    playState
}

mainCall <- function(playState = playDevCur()) {
    recursiveIndex(playState$call, playState$tmp$main.call.index)
}

"mainCall<-" <- function(playState = playDevCur(), value) {
    recursiveIndex(playState$call, playState$tmp$main.call.index) <- value
    playState
}

## used only by mainCall
recursiveIndex <- function(call, index) {
    ## if index is simple...
    if (length(index) == 1) {
        ## if index is NA, use original call object
        if (is.na(index)) return(call)
        return(call[[index]])
    }
    getx <- paste("[[", index, "]]", sep="", collapse="")
    eval(parse(text=paste("call", getx, sep="")))
}

## used only by "mainCall<-"
"recursiveIndex<-" <- function(call, index, value) {
    ## if index is simple...
    if (length(index) == 1) {
        ## if index is NA, use original call object
        if (is.na(index)) return(value)
        call[[index]] <- value
        return(call)
    }
    getx <- paste("[[", index, "]]", sep="", collapse="")
    eval(parse(text=paste("call", getx, " <- value", sep="")))
    call
}

updateMainCall <- function(playState = playDevCur()) {
    ## sets tmp$main.call.index, accepts.arguments, callName
    ## find which component of the call takes arguments (xlim etc)
    main.function <- playState$main.function
    tmpCall <- playState$call
    okCallPath <- function(tmpCall, main.function) {
        name <- toString(tmpCall[[1]])
        ## ignore expression() constructs (typically plotmath)
        #if (identical(name, quote(expression)))
        if (name %in% c("expression", "quote", "bquote", "substitute", "alist"))
            return(NULL)
        if (!is.null(main.function)) {
            ok <- identical(name, main.function)
        } else {
            if (is.symbol(tmpCall[[1]])) {
                tmpFun <- get(as.character(tmpCall[[1]]), mode = "function")
            } else {
                tmpFun <- eval(tmpCall[[1]])
            }
            ok <- any(c("xlim", "...") %in% names(formals(tmpFun)))
            ok <- ok && !(name == "with") ## skip `with` function
        }
        if (ok) return(TRUE)
        if (length(tmpCall) > 1)
            for (i in seq(2, length(tmpCall)))
                if (is.call(tmpCall[[i]])) {
                    tmpPath <- okCallPath(tmpCall[[i]], main.function)
                    if (isTRUE(tmpPath)) return(i)
                    if (!is.null(tmpPath)) return(c(i, tmpPath))
                }
        return(NULL)
    }
    main.call.index <-
        okCallPath(tmpCall, main.function)
    if (is.null(main.function)) {
        ## look for a call to "plot" ## TODO -- can drop this?
        main.call.index.plot <- okCallPath(tmpCall, "plot")
        if (!is.null(main.call.index.plot)) {
            ## found "plot" call
            main.call.index <- main.call.index.plot
        }
    }
    if (isTRUE(main.call.index)) main.call.index <- NA ## top-level
    ## check whether the called function accepts arguments
    playState$accepts.arguments <- !is.null(main.call.index)
    ## set index to top-level even if looks invalid, so callArg() works
    if (is.null(main.call.index)) main.call.index <- NA ## top-level
    playState$tmp$main.call.index <- main.call.index
    mainCall <- mainCall(playState)
    playState$callName <- toString(deparse(mainCall[[1]]))
    ## put call into canonical form, but with first argument un-named
    if (playState$accepts.arguments) {
        ## apply match.call()
        callFun <- eval(mainCall[[1]])
        firstArgName <- names(mainCall)[2]
        mainCall <- match.call(callFun, mainCall)
        if (is.null(firstArgName) || (firstArgName == ""))
            if (!is.null(names(mainCall))) names(mainCall)[2] <- ""
        mainCall(playState) <- mainCall
    }
}

rawXLim <- function(playState = playDevCur(), space="plot")
    rawXYLim(playState, space=space)$x

rawYLim <- function(playState = playDevCur(), space="plot")
    rawXYLim(playState, space=space)$y

rawXYLim <- function(playState, space="plot")
{
    playDevSet(playState)
    if (playState$is.lattice && (space == "plot")) {
        ## if space does not specify a panel, just pick one
        space <- packet.number()
        if (length(space) == 0) {
            packets <- playState$tmp$currentLayout
            space <- packets[packets > 0][1]
        }
        space <- paste("packet", space)
    }
    playDo(playState,
           list(x=convertX(unit(0:1, "npc"), "native", valueOnly=TRUE),
                y=convertY(unit(0:1, "npc"), "native", valueOnly=TRUE)),
           space=space)
}

"rawXLim<-" <- function(playState = playDevCur(), value)
{
    setRawXYLim(playState, value, "x")
    playState
}

"rawYLim<-" <- function(playState = playDevCur(), value)
{
    setRawXYLim(playState, value, "y")
    playState
}

setRawXYLim <- function(playState, x, x.or.y=c("x", "y"))
{
    playDevSet(playState)
    x.or.y <- match.arg(x.or.y)
     if (playState$is.lattice) {
         makeScalesArgAList(playState)
        ## TODO: packet 1 may not exist?
        x.panel <- xyData(playState, space="packet 1")[[x.or.y]]
        ## set factor labels explicitly, otherwise they are coerced to numeric
        if (is.factor(x.panel)) {
            scales.labels <- substitute(scales$w$labels,
                                        list(w = as.symbol(x.or.y)))
            scales.at <- substitute(scales$w$at,
                                        list(w = as.symbol(x.or.y)))
            if (is.null(callArg(playState, scales.labels))) {
                callArg(playState, scales.labels) <- levels(x.panel)
                callArg(playState, scales.at) <- 1:nlevels(x.panel)
            }
        }
        else if (is.somesortoftime(x.panel)) {
          class(x) <- class(x.panel)
          if (inherits(x.panel, "Date"))
            x <- call("as.Date", format(x))
          if (inherits(x.panel, "POSIXct"))
            x <- call("as.POSIXct", format(x))
          if (inherits(x.panel, "yearmon"))
            x <- call("as.yearmon", format(as.Date(x)))
          if (inherits(x.panel, "yearqtr"))
            x <- call("as.yearqtr", format(as.Date(x)))
        }
        else {
          ## numeric
          ## it seems now (lattice 0.17-12) that xlim/ylim are used directly
          ## so we don't need this:
          #isExtended <- switch(x.or.y,
          #                     x = playState$trellis$x.scales$axs == "r",
          #                     y = playState$trellis$y.scales$axs == "r")
          #f <- lattice.getOption("axis.padding")$numeric
          #if (isExtended) x <- shrinkrange(x, f=f)
        }
    }
    else if (!is.null(playState$viewport)) {
      ## non-lattice grid plot
      ## (do not know if the range is extended or not).
    }
    else {
      ## base graphics plot
      isExtended <- switch(x.or.y,
                           x = (par("xaxs") == "r"),
                           y = (par("yaxs") == "r"))
      if (isExtended) x <- shrinkrange(x, f=0.04)
    }
    ## convert back from log scale if required
    x <- spaceCoordsToDataCoordsXY(playState, x, x.or.y=x.or.y)
    ## round such that approximation error is within 1/1000 of x/y range
    if (is.numeric(x)) {
        digits <- max(3 - floor(log10(abs(diff(x)))), 0)
        x <- round(x, digits = digits)
    }
    if (x.or.y == "x") callArg(playState, "xlim") <- x
    if (x.or.y == "y") callArg(playState, "ylim") <- x
}

makeScalesArgAList <- function(playState)
{
    if (!isTRUE(playState$is.lattice)) return()
    scales <- callArg(playState, "scales")
    if (is.null(scales)) return()
    if (is.character(scales)) {
        callArg(playState, "scales") <-
            list(relation = scales)
        return()
    }
    if (is.character(scales$x)) {
        callArg(playState, quote(scales$x)) <-
            list(relation = scales$x)
    }
    if (is.character(scales$y)) {
        callArg(playState, quote(scales$y)) <-
            list(relation = scales$y)
    }
}

playSourceCode <- function(playState = playDevCur())
{
    theHeader <-
        paste("library(grid)",
              "library(lattice)",
              "library(playwith) ## (for panel.usertext, etc)",
              "## + might need others, often library(latticeExtra).",
              "## Assuming that the data are attached and any",
              "## customised style settings are in place; save with",
              "## myStyle <- trellis.par.get(); then restore with",
              "## trellis.par.set(myStyle)",
              sep = "\n")
    code <- list()
    comm <- list()
    code$plot <- playState$call
    if (playState$is.lattice) {
        code$plot <- call("print", playState$call)
        if (playState$pages > 1) {
            code$plot <- call("plotOnePage", playState$call,
                            page = playState$page)
        }
        ## use trellis$par.settings (if any) for annotations
        pars <- playState$trellis$par.settings
        if (length(pars) > 0) {
            code$plot <- c(code$plot,
                           call("<-", quote(opar),
                                call("trellis.par.set", pars)),
                           quote(on.exit(trellis.par.set(opar))))
        }

    }
    code$plot <- as.expression(code$plot)
    ## set up viewports
    comm$vps <- "set up viewports"
    code$vps <- expression()
    if (playState$is.lattice)
        code$vps <-
            c(code$vps,
              expression(
                  downViewport(trellis.vpname("toplevel"))
                  )
              )
    code$vps <-
        c(code$vps,
          expression(
              pushViewport(viewport(name = "pageAnnotationVp",
                                    yscale = c(1, 0))),
              upViewport(0))
          )
    if (playState$is.base) {
        code$vps <- c(code$vps, expression(
        local({
            vps <- baseViewports()
            vps$plot$name <- "plot"
            vps$plot$clip <- TRUE
            vps$plot.clip.off <-
                viewport(xscale=par("usr")[1:2],
                         yscale=par("usr")[3:4],
                         clip="off", name = "plot.clip.off")
            pushViewport(do.call("vpStack", vps))
        })))
    }
    ## annotations etc
    comm$linked <- "draw brushed (highlighted) points"
    code$linked <- drawLinkedLocal(playState, return.code = TRUE)
    comm$labels <- "add labels to data points"
    code$labels <- drawLabels(playState, return.code = TRUE)
    comm$annots <- "draw custom annotations"
    code$annots <- drawAnnotations(playState, return.code = TRUE)
    hasExtras <- with(code, (length(linked) || length(labels) ||
                             length(annots)))
    ## viewports are not needed unless drawing extras
    if (hasExtras == FALSE)
        code$vps <- NULL
    ## convert to text with interspersed comments
    theSource <- theHeader
    opts <- playwith.getOption("deparse.options")
    for (x in names(code)) {
        if (length(code[[x]]) > 0) {
            if (!is.null(comm[[x]]))
                theSource <- c(theSource,
                               paste("##", comm[[x]]))
            theSource <- c(theSource,
                           unlist(lapply(code[[x]], deparse, width = 42,
                                         control = opts)))
        }
    }
    ## clean up
    if (hasExtras)
        theSource <- c(theSource, "upViewport(0)")
    theSource <- paste(theSource, sep = "\n", collapse = "\n")
    theSource
}

playPrompt <- function(playState, text = NULL)
{
    if (is.null(text)) {
        playThawGUI(playState)
        playState$widgets$statusbar$pop(0)
    } else {
        playFreezeGUI(playState)
        playState$widgets$statusbar$push(0, toString(text))
    }
    invisible()
}

playFreezeGUI <- function(playState = playDevCur())
    playSetFreezeGUI(playState, TRUE)

playThawGUI <- function(playState = playDevCur())
    playSetFreezeGUI(playState, FALSE)

playSetFreezeGUI <- function(playState, frozen)
{
    playState$tmp$now.interacting <- frozen
    with(playState$widgets, {
        ## TODO: freeze GlobalActions etc?
        playState$actionGroups[["PlotActions"]]$setSensitive(!frozen)
        #topToolbar["sensitive"] <- !frozen
        #leftToolbar["sensitive"] <- !frozen
        #rightToolbar["sensitive"] <- !frozen
        ## leave bottom toolbar alone as this is where parameter
        ## control tools go (otherwise long slider drags interrupted).
        ## these tools check plot.ready before redrawing (threads...).
        ## similarly, leave page and time scrollbars as sensitive.
        if (!is.null(playState$widgets$latticist))
            latticist["sensitive"] <- !frozen
    })
    playState$win$getWindow()$setCursor(
      if (frozen) gdkCursorNew(GdkCursorType["watch"]) else NULL)
}

blockRedraws <- function(expr, playState = playDevCur())
{
    oval <- playState$tmp$skip.redraws
    playState$tmp$skip.redraws <- TRUE
    da <- playState$widgets$drawingArea
    daAlloc <- da$getAllocation()$allocation
    da$setSizeRequest(daAlloc$width, daAlloc$height)
                                        #playState$win$setGeometryHints(da, list(max.width=myW, min.width=myW,
                                        #	max.height=myH, min.height=myH))
                                        #da$window$freezeUpdates() # hmm
    foo <- try(eval.parent(substitute(expr)))
    ## try to force redraw
    gdkWindowProcessAllUpdates()
    while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)

                                        #da$window$thawUpdates()
                                        #playState$win$setGeometryHints(da, list())
    da$setSizeRequest(-1, -1)
    playState$tmp$skip.redraws <- oval
    foo
}

hideWidgetNoRedraw <- function(playState, widget, horiz)
{
    whichDim <- if (horiz) "height" else "width"
    if (widget["visible"]) blockRedraws({
        widgSize <- widget$getAllocation()$allocation
        winSize <- playState$win$getSize()
        widget["visible"] <- FALSE
        winSize[[whichDim]] <- winSize[[whichDim]] - widgSize[[whichDim]]
        playState$win$resize(winSize$width, winSize$height)
    })
}

## TODO: store value in playState
isBasicDeviceMode <- function(playState)
{
    if ((length(playState$call) == 1) &&
        identical(playState$call[[1]], quote(`{`))) {
        ## basic device mode
        ## (do not know the call)
        return(TRUE)
    }
    FALSE
}

