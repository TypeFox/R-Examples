## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

globalActionGroup <- function(playState)
{
    entries <-
        list( ## : name, stock icon, label, accelerator, tooltip, callback
             list("Clone", "gtk-new", "_New linked plot", "<Ctrl>N", "Make a linked clone of this window (with sync'd brushing)", clone_handler),
             list("Unlink", "gtk-disconnect", "_Unlink", NULL, "Remove links to other plots (sync'd brushing)", unlink_handler),
             list("Save", "gtk-save", "_Save", "<Ctrl>S", "Export current plot to an image file", save_handler),
             list("Copy", "gtk-copy", "_Copy", "<Ctrl><Shift>C", "Copy current plot as an image", copy_handler),
             list("Print", "gtk-print", "_Print", "<Ctrl>P", "Print current plot", print_handler),
             list("Close", "gtk-close", "Close", "<Ctrl>W", "Close window and device", close_handler),
             list("SetSize", NULL, "Set device _size...", "<Ctrl>M", NULL, set.size_handler),
             list("IncrFont", NULL, "_Increase font size", "<Ctrl>plus", NULL, incr.font_handler),
             list("DecrFont", NULL, "De_crease font size", "<Ctrl>minus", NULL, decr.font_handler),
             list("StyleSettings", "gtk-select-color", "_Style settings", "<Ctrl>slash", NULL, style.settings_handler),
             list("EditCall", "gtk-edit", "_Edit call...", "<Ctrl>E", "Edit the plot call", edit.call_handler),
             list("Back", "gtk-go-back", "Back", "<Alt>Left", "Go back to previous plot call", back_handler),
             list("Forward", "gtk-go-forward", "Forward", "<Alt>Right", "Go to next plot call", forward_handler),
             list("Redraw", "gtk-refresh", "Re_draw", "<Ctrl>R", NULL, redraw_handler),
             list("Reload", "gtk-refresh", "_Reload and redraw", "<Ctrl><Shift>R", NULL, reload_handler),
             list("Interrupt", "gtk-stop", "_Stop", "<Ctrl>period", "Stop (interrupt) a plot", interrupt_handler),
             list("SaveCode", "gtk-save", "Save c_ode", "<Ctrl><Shift>S", "Save R code for this plot", save.code_handler),
             list("ViewSource", NULL, "Plot s_ource", "<Ctrl>U", NULL, view.source_handler),
             list("HelpPlot", "gtk-help", "_Help for this plot", "F1", "Open help page for this plot command", help_handler),
             list("HelpPlaywith", NULL, "help(playwith)", NULL, NULL, help.playwith_handler),
             list("About", NULL, "_About playwith", NULL, NULL, about_handler),
             list("Website", NULL, "_Website", NULL, "The playwith website (for contact, bugs, etc)", website_handler),
             list("SummariseData", NULL, "Summarise object"),
             list("EditData", NULL, "Edit object"),
             list("SaveData", NULL, "Save object"),
             list("PurgeObjects", NULL, "Purge object cache")
             )

    toggleEntries <-
        list( ## : name, stock icon, label, accelerator, tooltip, callback, active?
             list("Keep", "gtk-media-stop", "_Keep open", "<Ctrl>D", "Keep this plot (do not replace with the next plot)", keep_handler, FALSE),
             list("StayOnTop", "gtk-leave-fullscreen", "St_ay on top", "<Ctrl>grave", "Show this window above all others", stay.on.top_handler, FALSE),
             ## options (uiOptionsActions.R)
             list("ClipAnnot", NULL, "_Clip annotations", NULL, "", clip.annotations_handler, FALSE),
             list("PageAnnot", NULL, "_Annot. on page (fixed pos.)", NULL, "Place annotations with respect to the page, not plot coordinates", page.annotation_handler, FALSE),
             list("ShowStatusbar", NULL, "Status _bar", NULL, NULL, show.statusbar_handler, TRUE),
             list("ShowToolbars", NULL, "Toolbars", NULL, NULL, show.toolbars_handler, TRUE),
             list("ShowTooltips", NULL, "T_ooltips", NULL, "", show.tooltips_handler, FALSE)
             )

    ## construct action group with playState passed to callbacks
    aGroup <- gtkActionGroupNew("GlobalActions")
    aGroup$addActions(entries, playState)
    aGroup$addToggleActions(toggleEntries, playState)
    aGroup
}

updateGlobalActions <- function(playState)
{
    aGroup <- playState$actionGroups[["GlobalActions"]]
    ## Unlink
    hasLinks <- (length(playState$linked$subscribers) > 1)
    aGroup$getAction("Unlink")$setVisible(hasLinks)
    ## Back, Forward
    callEntry <- playState$widgets$callEntry
    nHistory <- callEntry$getModel()$iterNChildren()
    aGroup$getAction("Forward")$setSensitive(callEntry["active"] > 0)
    aGroup$getAction("Back")$setSensitive(callEntry["active"] < nHistory-1)
    ## Keep
    aGroup$getAction("Keep")$setActive(isTRUE(playState$keep))
    ## StayOnTop
    aGroup$getAction("StayOnTop")$setActive(isTRUE(playState$stay.on.top))
}

clone_handler <- function(widget, playState)
{
    playDevSet(playState)
    ## start off with an empty shell
    sizein <- dev.size("in")
    newOne <- playwith({}, new = TRUE,
                       title = paste("Clone of", playState$win["title"]),
                       width = sizein[1], height = sizein[2],
                       pointsize = playState$pointsize)
    ## copy everything except widgets / pointers / ID / tmp
    for (name in ls(playState)) {
        if (name %in% c("win", "dev", "tmp", "ID",
                        "uiManager", "actionGroups", "widgets"))
            next
        ## copy to the clone
        assign(name, get(name, playState), newOne)
    }
    ## link it
    i <- length(playState$linked$subscribers) + 1
    playState$linked$subscribers[[i]] <- newOne
    playNewPlot(newOne)
    updateGlobalActions(playState)
}

unlink_handler <- function(widget, playState)
    playUnlink(playState)

save_handler <- function(widget, playState)
{
    ## disable toolbars until this is over
    playFreezeGUI(playState)
    on.exit(playThawGUI(playState))
    ## get filename
    myExt <- playwith.getOption("save.as.format")
    myDefault <- if (!is.null(playState$title))
        playState$title else paste(playState$callName, "%02d", sep="")
    myDefault <- paste(myDefault, myExt, sep=".")
    ## construct save file dialog
    okExt <- c("pdf","png","jpg","jpeg","tif","tiff","ps","eps","svg","wmf","emf","fig")
    filter <- list("All files" = list(patterns = c("*")),
                   "PDF" = list(patterns = c("*.pdf")),
                   "PNG (bitmap)" = list(patterns = c("*.png")),
                   "EPS" = list(patterns = c("*.eps")),
                   "SVG" = list(patterns = c("*.svg")),
                   "WMF (metafile)" = list(patterns = c("*.wmf", "*.emf")),
                   "JPEG" = list(patterns = c("*.jpg", "*.jpeg")),
                   "TIFF" = list(patterns = c("*.tif", "*.tiff")),
                   "xfig" = list(patterns = c("*.fig")))
    ## TODO: pdfWriter with bitmap()
    filename <- gfile("Export plot to image file", type = "save",
                      filter = filter, initialfilename = myDefault)
    if (is.na(filename)) return()
    ext <- tolower(get.extension(filename))
    if ((ext %in% okExt) == FALSE) {
        filename <- paste(filename, myExt, sep=".")
        ext <- myExt
    }
    ## save plot to file
    playDevSet(playState)
    ## note: baseViewports will be corrupted if device size changes
    ## so need to keep the same size with dev.copy()...
    w.in <- dev.size("in")[1]
    h.in <- dev.size("in")[2]
    ## use same pointsize as embedded device
    ps <- playState$pointsize
    #if (ext == "eps") setEPS()
    devName <-
        switch(ext,
               pdf = "pdf",
               ps =, eps = "postscript",
               png = "Cairo_png",
               jpeg =, jpg = "jpeg",
               tiff =, tif = "tiff",
               svg = "Cairo_svg",
               wmf =, emf = "win.metafile",
               fig = "xfig",
           {gmessage.error("Unrecognised filename extension")
            stop("Unrecognised filename extension")})
    devFun <- get(devName, mode = "function")
    callNm <- if (ext == "eps") "dev.copy2eps" else "dev.copy"
    devCall <- call(callNm, device = devFun, file = filename,
                    width = w.in, height = h.in, pointsize = ps)
    if (!is.null(formals(devFun)$units))
        devCall$units <- "in"
    eval(devCall)
    dev.off()
}

copy_handler <- function(widget, playState)
{
    ## disable toolbars until this is over
    playFreezeGUI(playState)
    on.exit(playThawGUI(playState))
    playDevSet(playState)
    w.in <- dev.size("in")[1]
    h.in <- dev.size("in")[2]
    ## use same pointsize as embedded device
    ps <- playState$pointsize
    if (exists("win.metafile")) { ## i.e. in MS Windows
        copy.exts <- c("wmf", "png")
        copy.labels <- c("Windows Metafile (wmf)", "Bitmap (png)")
        sel.label <- select.list(copy.labels, preselect=copy.labels[1],
                                 title="Copy in which format?")
        playState$win$present()
        if (sel.label == "") return()
        copy.ext <- copy.exts[copy.labels == sel.label]
        if (copy.ext == "wmf") {
            dev.copy(win.metafile, width=w.in, height=h.in,
                     pointsize = ps)
            dev.off()
            playState$win$present()
            return()
        }
    }
    ## save plot to temporary file, to copy as png
    filename <- paste(tempfile(), ".png", sep="")
    dev.copy(Cairo_png, file=filename, width=w.in, height=h.in,
             pointsize = ps)
    dev.off()
    im <- gdkPixbufNewFromFile(filename)$retval
    gtkClipboardGet("CLIPBOARD")$setImage(im)
    file.remove(filename)
}

print_handler <- function(widget, playState)
{
    playDevSet(playState)
    isWindows <- (.Platform$OS.type == "windows")
    ## note dev.print copies the current width / height / ps
    if (isWindows) dev.print(win.print)
    else dev.print()
}

close_handler <- function(widget, playState)
    window.close_handler(playState = playState)

edit.call_handler <- function(widget, playState)
{
    callTxt <-
        paste(deparse(playState$call, width.cutoff = 42, control=
                      playwith.getOption("deparse.options")),
              collapse="\n")
    repeat {
        newTxt <- guiTextInput(callTxt, title="Edit plot call",
                               prompt="", accepts.tab = FALSE)
        if (is.null(newTxt)) break
        callTxt <- newTxt
        tmp <- tryCatch(parse(text=callTxt), error=function(e)e)
        ## check whether there was a syntax error
        if (inherits(tmp, "error")) {
            gmessage.error(conditionMessage(tmp))
        } else {
            ## if more than one call, wrap them in braces
            playState$call <- if (length(tmp) > 1)
                as.call(c(as.symbol("{"), tmp)) else tmp[[1]]
            playNewPlot(playState)
            break
        }
    }
    playState$win$present()
}

edit.call.inline_handler <- function(widget, playState)
{
    ## the original call
    callTxt <- deparseOneLine(playState$call, control=
                              playwith.getOption("deparse.options"))
    newTxt <- widget["text"]
    if (identical(newTxt, callTxt)) return()
    if (identical(newTxt, "")) return()
    tmp <- tryCatch(parse(text=newTxt), error=function(e)e)
    ## check whether there was a syntax error
    if (inherits(tmp, "error")) {
        gmessage.error(conditionMessage(tmp))
    } else {
        ## if more than one call, wrap them in braces
        playState$call <- if (length(tmp) == 1) tmp[[1]]
            else as.call(c(as.symbol("{"), tmp))
        playNewPlot(playState)
    }
    playState$win$present()
}

set.size_handler <- function(widget, playState) {
    da <- playState$widgets$drawingArea
    daAlloc <- da$getAllocation()$allocation
    owidth <- daAlloc$width
    oheight <- daAlloc$height
    ## prompt user for new size
    widthW <- gedit(toString(owidth), width=7, coerce.with=as.numeric)
    heightW <- gedit(toString(oheight), width=7, coerce.with=as.numeric)
    unitVals <- c("pixels", "cm", "inches", "percent")
    ## functions to convert between selected units and pixels
    dpi <- dev.size("px")[1] / dev.size("in")[1]
    px2a <- function(px, unit, ref.px)
        signif(switch(unit, pixels = round(px), percent = 100 * (px/ref.px),
                      inches = px/dpi, cm = 2.54 * px/dpi), 3)
    a2px <- function(a, unit, ref.px)
        round(switch(unit, pixels = a, percent = ref.px * (a/100),
                      inches = a*dpi, cm = (a / 2.54)*dpi))
    unitsW <- gcombobox(unitVals, handler = function(h, ...) {
        wnum <- px2a(owidth, svalue(h$obj), owidth)
        hnum <- px2a(oheight, svalue(h$obj), oheight)
        svalue(widthW) <- toString(wnum)
        svalue(heightW) <- toString(hnum)
        })
    lay <- glayout()
    lay[1,1] <- "Width: "
    lay[2,1] <- "Height: "
    lay[1,2] <- widthW
    lay[2,2] <- heightW
    lay[2,3] <- unitsW
    width <- height <- NA
    result <- gbasicdialog(title="Set device size", widget=lay,
                           handler=function(...) {
                               unit <- svalue(unitsW)
                               width <<- a2px(svalue(widthW), unit, owidth)
                               height <<- a2px(svalue(heightW), unit, oheight)
                           })
    playState$win$present()
    if (!isTRUE(result)) return()
    width <- max(10, width)
    height <- max(10, height)
    da$setSizeRequest(width, height)
    playState$win$resize(1, 1) ## as small as possible
    ## try to force resize
    gdkWindowProcessAllUpdates()
    while (gtkEventsPending()) gtkMainIterationDo(blocking=FALSE)
    ## remove constraint after resize
    da$setSizeRequest(-1, -1)
}

incr.font_handler <- function(widget, playState)
{
    playDevSet(playState)
    ps <- playState$pointsize <- playState$pointsize + 1
    par(ps = ps)
    trellis.par.set(fontsize = list(text = ps))
    playReplot(playState)
}

decr.font_handler <- function(widget, playState)
{
    playDevSet(playState)
    ps <- playState$pointsize <- playState$pointsize - 1
    par(ps = ps)
    trellis.par.set(fontsize = list(text = ps))
    playReplot(playState)

}

style.settings_handler <- function(widget, playState) {
    if (!require("latticist")) {
        gmessage(paste('This feature requires the latticist package.',
                       'Install it with install.packages("latticist")'))
        return()
    }
    playDevSet(playState)
    isBase <- (!playState$is.lattice && is.null(playState$viewport))
    latticeStyleGUI(pointsize = playState$pointsize,
                    base.graphics = isBase)
}

back_handler <- function(widget, playState) {
    with(playState$widgets,
        callEntry["active"] <- callEntry["active"] + 1)
}

forward_handler <- function(widget, playState) {
    with(playState$widgets,
        callEntry["active"] <- callEntry["active"] - 1)
}

redraw_handler <- function(widget, playState)
    playReplot(playState)

reload_handler <- function(widget, playState)
    playNewPlot(playState)

interrupt_handler <- function(widget, playState) {
    .C(do_interrupt)
}

save.code_handler <- function(widget, playState)
{
    myDefault <- "plot.R"
    filter <- list("All files" = list(patterns = c("*")),
                   "Text files" = list(patterns = c("*.txt")),
                   "R files" = list(patterns = c("*.R")))
    filename <- gfile("Export plot source code", type = "save",
                      filter = filter, initialfilename = myDefault)
    if (is.na(filename)) return()
    theSource <- playSourceCode(playState)
    cat(theSource, file = filename)
}

view.source_handler <- function(widget, playState)
{
    theSource <- playSourceCode(playState)
    guiTextView(theSource, title = "Plot source code")
}

help_handler <- function(widget, playState)
{
    if (playState$accepts.arguments == FALSE) {
        gmessage.error("Do not know the name of the plot function.")
        return()
    }
    ## work out which (S3) method was called, if any
    callName <- playState$callName
    methNames <- methods(callName)
    if (length(methNames) > 0) {
        myClass <- try(class(callArg(playState, 1)), silent=TRUE)
        if (!inherits(myClass, "try-error")) {
            myMeth <- paste(callName, myClass, sep=".")
            ok <- (myMeth %in% methNames)
            if (any(ok)) callName <- myMeth[ok][1]
        }
    }
    print(help(callName))
}

help.playwith_handler <- function(widget, playState)
    print(help("playwith"))

about_handler <- function(widget, playState) {
    activate.email <- function(about, link, data)
        browseURL(paste("mailto:", link, sep=""))
    activate.url <- function(about, link, data)
        browseURL(link)
    gtkAboutDialogSetEmailHook(activate.email)
    gtkAboutDialogSetUrlHook(activate.url)

    ## TODO: this shows the wrong name! ("Rterm.exe")
    gtkShowAboutDialog(playState$win,
                       title = "playwith",
                       name = "playwith",
                       version = packageDescription("playwith")$Version,
                       comments = "interactive plots in R using GTK+",
                       copyright = "(C) 2008 Felix Andrews",
                       license = "GNU General Public Licence version 2 or later",
                       website = "http://playwith.googlecode.com/",
                       authors = c("Felix Andrews <felix@nfrac.org>"),
                       documenters = c("Felix Andrews")
                       )
}

website_handler <- function(...)
    browseURL("http://playwith.googlecode.com/")

keep_handler <- function(widget, playState)
    playState$keep <- widget["active"]

stay.on.top_handler <- function(widget, playState) {
    playState$stay.on.top <- widget["active"]
    playState$win$setKeepAbove(widget["active"])
}

