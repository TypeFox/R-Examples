## playwith: interactive plots in R using GTK+
##
## Copyright (c) 2007 Felix Andrews <felix@nfrac.org>
## GPL version 2 or newer

parameterControlTool <-
    function(playState, name, value,
             label = name, handler = NULL,
             horizontal = TRUE)
{
    stopifnot(length(value) > 0)
    if (is.list(value)) {
        if (!is.null(value$label))
            label <- value$label
        if (!is.null(value$handler))
            handler <- value$handler
        value <- value[[1]]
    }
    if (is.character(handler))
        handler <- get(handler)
    if (!is.logical(value) && !is.function(value))
        label <- paste(label, ": ", sep="")
    ## signal handlers
    updateParamValue <- function(widget, playState) {
        ## (can't use widget["value"] with a slider)
        newval <- widget$getValue()
        oldval <- get(name, envir=playState$env)
        if (identical(oldval, newval)) return()
        assign(name, newval, envir=playState$env)
        ## run custom 'handler', abort if returns FALSE
        if (!is.null(handler)) {
            result <- handler(playState, newval)
            if (identical(result, FALSE)) return()
        }
        if (!isTRUE(playState$tmp$plot.ready)) return()
        playReplot(playState)
    }
    updateParamText <- function(widget, playState) {
        newval <- widget["text"]
        oldval <- get(name, envir=playState$env)
        if (identical(oldval, newval)) return()
        assign(name, newval, envir=playState$env)
        ## run custom 'handler', abort if returns FALSE
        if (!is.null(handler)) {
            result <- handler(playState, newval)
            if (identical(result, FALSE)) return()
        }
        if (!isTRUE(playState$tmp$plot.ready)) return()
        playReplot(playState)
    }
    updateParamTextNumeric <- function(widget, playState) {
        newval <- try(as.numeric(widget["text"]))
        if (inherits(newval, "try-error")) return()
        if (is.na(newval)) return()
        oldval <- get(name, envir=playState$env)
        if (identical(oldval, newval)) return()
        assign(name, newval, envir=playState$env)
        ## run custom 'handler', abort if returns FALSE
        if (!is.null(handler)) {
            result <- handler(playState, newval)
            if (identical(result, FALSE)) return()
        }
        if (!isTRUE(playState$tmp$plot.ready)) return()
        playReplot(playState)
    }
    updateParamCombobox <- function(widget, playState) {
        ## signal also emitted on typing, ignore
        if (widget$getActive() == -1) return()
        newval <- widget$getActiveText()
        oldval <- get(name, envir=playState$env)
        if (identical(oldval, newval)) return()
        assign(name, newval, envir=playState$env)
        ## run custom 'handler', abort if returns FALSE
        if (!is.null(handler)) {
            result <- handler(playState, newval)
            if (identical(result, FALSE)) return()
        }
        if (!isTRUE(playState$tmp$plot.ready)) return()
        playReplot(playState)
    }
    updateParamActive <- function(widget, playState) {
        newval <- widget["active"]
        oldval <- get(name, envir=playState$env)
        if (identical(oldval, newval)) return()
        assign(name, newval, envir=playState$env)
        ## run custom 'handler', abort if returns FALSE
        if (!is.null(handler)) {
            result <- handler(playState, newval)
            if (identical(result, FALSE)) return()
        }
        if (!isTRUE(playState$tmp$plot.ready)) return()
        playReplot(playState)
    }
    ## construct widget based on the type of value.
    ## note that the initial value has been set in playwith()
    ## integer or AsIs : spinbutton
    if (is.integer(value) || inherits(value, "AsIs")) {
        if (inherits(value, "AsIs")) value <- as.vector(value)
        box <- gtkVBox()
        box$packStart(gtkLabel(label))
        if (length(value) == 1) {
            ## only one value given -- make up a range
            range <- 10 * (1 + abs(value)) ^ 2 * c(-1, 1)
        } else {
            ## range given
            range <- range(value)
        }
        step <- min(diff(unique(sort(value))))
        widget <- gtkSpinButton(min=min(range), max=max(range), step=step)
        widget["digits"] <- max(0, - floor(log10(step)))
        widget$setValue(get(name, envir=playState$env))
        gSignalConnect(widget, "value-changed",
                       updateParamValue, data=playState)
        box$packStart(widget)
        foo <- gtkToolItem()
        foo$add(box)
        return(foo)
    }
    ## numeric: entry or slider
    if (is.numeric(value)) {
        if (length(value) == 1) {
            ## entry coercing to numeric
            box <- gtkVBox()
            box$packStart(gtkLabel(label))
            widget <- gtkEntry()
            widget["text"] <- toString(get(name, envir=playState$env))
            widget["width-chars"] <- 6
            gSignalConnect(widget, "activate",
                           updateParamTextNumeric, data=playState)
            box$packStart(widget)
            foo <- gtkToolItem()
            foo$add(box)
            return(foo)
            ## only one value given -- make up a range
#            range <- 10 * (1 + abs(value)) ^ 2 * c(-1, 1)
#            step <- 10^round(log10(abs(value))-1)
#            if (value == 0) step <- 1
        }
        ## scale (slider)
        if (length(value) == 2) {
            ## range given -- make up a step
            range <- range(value)
            step <- 10^round(log10(diff(range))-1)
        } else {
            ## explicit sequence given (but might not be regular)
            range <- range(value)
            step <- median(abs(diff(sort(value))))
        }
        digits <- max(0, - floor(log10(step)))
        box <- if (horizontal) gtkHBox() else gtkVBox()
        box$packStart(gtkLabel(label), expand=FALSE)
        widget <- if (horizontal)
            gtkHScale(min=min(range), max=max(range), step=step)
        else gtkVScale(min=min(range), max=max(range), step=step)
        widget$setValue(get(name, envir=playState$env))
        widget["digits"] <- digits
        widget["update-policy"] <- GtkUpdateType["delayed"]
        gSignalConnect(widget, "value-changed",
                       updateParamValue, data=playState)
        box$packStart(widget)
        foo <- gtkToolItem()
        foo$setExpand(TRUE)
        foo$add(box)
        return(foo)
    }
    ## character: entry or combobox
    if (is.character(value)) {
        box <- gtkVBox()
        box$packStart(gtkLabel(label))
        if (length(value) == 1) {
            ## text entry
            widget <- gtkEntry()
            widget["text"] <- get(name, envir=playState$env)
            #widget["width-chars"] <- 30
            gSignalConnect(widget, "activate",
                           updateParamText, data=playState)
        } else {
            ## combo box
            widget <- gtkComboBoxEntryNewText()
            widget$show()
            for (item in value)
                widget$appendText(item)
            #widget$setActiveText(get(name, envir=playState$env))
            index <- match(get(name, envir=playState$env), value)
            if (is.na(index)) index <- 1
            widget["active"] <- (index - 1)
            #gSignalConnect(widget, "editing-done",
            gSignalConnect(widget$getChild(), "activate",
                           updateParamText, data=playState)
            gSignalConnect(widget, "changed",
                           updateParamCombobox, data=playState)
        }
        box$packStart(widget)
        foo <- gtkToolItem()
        foo$add(box)
        return(foo)
    }
    ## logical: checkbutton
    if (is.logical(value)) {
        ## toggle button / checkbox
        widget <- gtkCheckButton(label)
        widget["active"] <- isTRUE(get(name, envir=playState$env))
        gSignalConnect(widget, "clicked",
                       updateParamActive, data=playState)
        foo <- gtkToolItem()
        foo$add(widget)
        return(foo)
    }
    ## function: button
    if (is.function(value)) {
        widget <- gtkButton(label)
        gSignalConnect(widget, "clicked",
                       function(widget, playState)
                       value(playState),
                       data = playState)
        foo <- gtkToolItem()
        foo$add(widget)
        return(foo)
    }
    ## otherwise...
    stop("do not know about ",
         toString(class(value)), " objects")
}
