#guiPane.tcltk <- function (dlg, item, ...)
#{
#    ## This is a basic pane
#
#    ## Check arguments
#    if (!inherits(dlg, "guiDlg"))
#        stop("'dlg' must be a 'guiDlg' object!")
#    if (!is.numeric(item) || item < 1 || item > length(dlg$panes))
#        stop("'item' is not numeric or is out of range!")
#    item <- as.integer(item)
#
#    ## Get general dialog parameters needed to organize the pane
#    cnt <- dlg$container
#    style <- dlg$style
#    width <- dlg$width
#    pixwidth <- dlg$pixwidth
#    labelwidth <- dlg$labelwidth
#    pane <- dlg$panes[[item]]
#
#    ## Get my pane parameters: argname, message & type
#    argname <- pane$argname
#    message <- pane$message
#    if (!is.null(message)) message <- paste(pane$message, collapse = "\n")
#    type <- pane$type
#    if (is.null(type)) type <- "entry"  # Default type
#
#    ## Construct the pane
#    ## Rem: button justify does not work => use fixed font and pad " " at the end!
#    paneFrame <- tkframe(cnt)
#    if (!is.null(argname)) {
#        ## Determine if this argument is included (by default, not!)
#        varUseIt <- tclVar("0")
#        arglab <- paste(argname, "=", sep = "")
#        argl <- nchar(arglab)
#        if (argl < labelwidth) {
#            arglab <- paste(arglab, paste(rep(" ", labelwidth - argl),
#            collapse = ""), sep = "")
#            argl <- labelwidth
#        }
#
#        butArg <- tkbutton(paneFrame, text = arglab,
#        font = style$font.fixed, fg = style$fg[5], width = argl,
#        takefocus = "0", relief = "flat", overrelief = "raised")
#        butArgToggle <- butArg
#        onArg <- function () {
#            ## Toggle varUseIt and color (emphasized <-> selected!)
#            if (tclvalue(varUseIt) == "0") {  # Select
#                tclvalue(varUseIt) <- "1"
#                tkconfigure(butArgToggle, fg = style$fg[4])
#            } else {  # Deselect
#                tclvalue(varUseIt) <- "0"
#                tkconfigure(butArgToggle, fg = style$fg[5])
#            }
#        }
#        tkconfigure(butArg, command = onArg)
#        ## This is used when the argument is edited
#        onArgEdit <- function () {
#            if (tclvalue(varUseIt) == "0") {
#                tclvalue(varUseIt) <- "1"
#                tkconfigure(butArgToggle, fg = style$fg[4])
#            }
#            return(tclVar(TRUE))  # Must return TRUE to accept edition!
#        }
#    }
#    if (!is.null(message)) {
#        labMessage <- tklabel(paneFrame, text = message, font = style$font.label,
#        justify = "left", fg = style$fg[1], wraplength = as.integer(
#        pixwidth / (labelwidth + width) * width))
#        if (!is.null(argname)) tkgrid(butArg, labMessage) else tkgrid(labMessage)
#        tkgrid(paneFrame, sticky = "w", padx = style$pads[1])
#        paneFrame <- tkframe(cnt)  # New frame
#        butArg <- tkbutton(paneFrame, text = " ", width = labelwidth,
#        font = style$font.fixed, relief = "flat", state = "disabled")
#    } else {  # There is no message... just display argname if not NULL
#        if (!is.null(argname)) {
#            if (argl > labelwidth) {  # arg too large, display it on top
#                tkgrid(butArg)
#                tkgrid(paneFrame, sticky = "w", padx = style$pads[1])
#                paneFrame <- tkframe(cnt)
#                butArg <- tkbutton(paneFrame, text = " ", width = labelwidth,
#                font = style$font.fixed, relief = "flat", state = "disabled")
#            }
#        } else {  # Neither message, nor label
#            if (labelwidth > 0)
#            butArg <- tkbutton(paneFrame, text = " ", width = labelwidth,
#            font = style$font.fixed, relief = "flat", state = "disabled")
#        }
#    }
#    if (labelwidth == 0) butArg <- NULL
#
#    ## Call guiPane.<type>.tcltk() functions to install specific widgets
#    fun <- paste("guiPane", type, "tcltk", sep = ".")
#    if (!exists(fun, where = 1, mode = "function")) fun <- "guiPane.entry.tcltk"
#    resenv <- get(fun, pos = 1, mode = "function")(
#    paneFrame, butArg, onArgEdit, varUseIt, dlg, item, ...)
#    if (is.null(resenv) && fun != "guiPane.entry.tcltk")  # Try default one
#    resenv <- get("guiPane.entry.tcltk", pos = 1, mode = "function")(
#    paneFrame, butArg, onArgEdit, varUseIt, dlg, item, ...)
#    ## Record resenv in dlg
#    dlg$panes[[item]]$env <- resenv
#    ## Place paneFrame in the dialog box
#    tkgrid(paneFrame, padx = style$pads[1], pady = style$pads[3], sticky = "w")
#
#    ## Return the modified dlg object
#    return(invisible(dlg))
#}
#
#guiPane.entry.tcltk <- function (paneFrame, butArg, onArgEdit, varUseIt, dlg,
#item, ...)
#{
#    ## This is a simple text entry pane
#
#    ## Get general dialog parameters needed here
#    style <- dlg$style
#    width <- dlg$width
#    pane <- dlg$panes[[item]]
#
#    ## Get my pane parameters: argname, default & fixedfont
#    argname <- pane$argname
#    default <- pane$default[1]
#    if (is.null(default) || is.na(default)) default <- "" else
#    default <- as.character(default)
#    fixedfont <- pane$fixedfont
#    if (is.null(fixedfont) || is.na(fixedfont)) fixedfont <- FALSE else
#    fixedfont <- (fixedfont == TRUE)
#
#    ## Install the specific widgets
#    varText <- tclVar(default)
#    if (fixedfont) Font <- style$font.fixed else Font <- style$font.text
#    txt <- tkentry(paneFrame, textvariable = varText, width = width,
#    font = Font, fg = style$fg[1], background = "white", relief = style$relief)
#    if (is.null(butArg)) tkgrid(txt) else tkgrid(butArg, txt)
#    if (!is.null(argname))
#        tkconfigure(txt, validate = "key", validatecommand = onArgEdit)
#    tkselection.from(txt, "0")
#    tkselection.to(txt, "end")
#    tkicursor(txt, "end")
#
#    ## Define the environment to interact with these widgets
#    resenv <- new.env(parent = parent.frame(2))
#    assign("varText", varText, envir = resenv)
#    assign("txt", txt, envir = resenv)
#    assign("argname", argname, envir = resenv)
#    result <- function() {
#    res <- if (is.null(argname)) tclvalue(varText) else
#    if (tclvalue(varUseIt) == "1") {
#        ## Compute the code for argument
#        if (argname == "...") tclvalue(varText) else
#        paste(argname, "=", tclvalue(varText))
#    } else ""
#        res
#    }
#    assign("result", result, envir = resenv)
#    select <- function() {
#        tkselection.from(txt, "0")
#        tkselection.to(txt, "end")
#        tkicursor(txt, "end")
#        tkfocus(txt)
#    }
#    assign("select", select, envir = resenv)
#    return(resenv)
#}
#
#guiPane.list.tcltk <- function (paneFrame, butArg, onArgEdit, varUseIt, dlg,
#item, ...)
#{
#    ## This is a single selection listbox
#
#    ## Get general dialog parameters needed here
#    style <- dlg$style
#    width <- dlg$width
#    pane <- dlg$panes[[item]]
#
#    ## Get my pane parameters: argname, choices, default, sort, listheight
#    argname <- pane$argname
#    choices <- pane$choices
#    if (is.null(choices)) choices = ""  # To make sure it works all the time
#    N <- length(choices)
#    if (!inherits(choices, "character") && N < 1)
#        stop("Pane", item, ": 'choices' must be a vector of strings!")
#    default <- pane$default[1]
#    if (!is.null(default)) {
#        if (!is.numeric(default))
#            stop("Pane", item, ": 'default' must be numeric or NULL!")
#        default <- as.integer(default)
#        if (default < 1 || default > N)
#            stop("Pane", item, ": 'default' is outside range!")
#    }
#    sort <- pane$sort
#    if (!is.null(sort) && !is.na(sort)) sort <- (sort == TRUE) else sort <- FALSE
#    if (sort) {  # Sort choices alphabetically
#        if (!is.null(default)) default <-  (1:N)[match(default, order(choices))]
#        choices <- sort(choices)
#    }
#    height <- pane$listheight
#    if (is.null(height)) height <- 4  # Default value
#    if (!is.numeric(height) || height < 2)
#        stop("Pane", item, ": 'listheight' must be numeric and > 1!")
#    height <- as.integer(height)
#
#    ## Install the specific widgets
#    scr <- tkscrollbar(paneFrame, repeatinterval = 5,
#    command = function(...) tkyview(tl, ...))
#    tl <- tklistbox(paneFrame, width = width - 2, height = height,
#    selectmode = "browse", yscrollcommand = function(...) tkset(scr, ...),
#    font = style$font.text, fg = style$fg[1], background = "white",
#    relief = style$relief, activestyle = "dotbox")
#    if (is.null(butArg)) tkgrid(tl, scr) else {
#        tkgrid(butArg, tl, scr)
#        tkgrid.configure(butArg, sticky = "nw")
#    }
#    tkgrid.configure(scr, rowspan = 5, sticky = "nsw")
#    tkgrid(paneFrame, padx = style$pads[1], pady = style$pads[3], sticky = "w")
#    for (i in 1:(length(choices)))
#        tkinsert(tl, "end", choices[i])
#    if (!is.null(default)) {
#        for (i in 1:length(default))
#            tkselection.set(tl, default[i] - 1)
#        tkyview(tl, default[1] - 1)  # Make sure selected item is visible
#    }
#    if (!is.null(argname))
#        tkbind(tl, "<<ListboxSelect>>", onArgEdit)
#
#    ## Define the environment to interact with these widgets
#    resenv <- new.env(parent = parent.frame(2))
#    assign("choices", choices, envir = resenv)
#    assign("tl", tl, envir = resenv)
#    assign("argname", argname, envir = resenv)
#    result <- function() {
#    sel <- choices[as.numeric(tkcurselection(tl)) + 1]
#    res <- if (is.null(argname)) paste(sel, collapse = ", ") else
#        if (tclvalue(varUseIt) == "1") {
#        ## Compute the code for argument
#        if (is.null(sel) || length(sel) == 0) sel <- "NULL"
#        if (length(sel) > 1) sel <- paste("c(", paste(sel, collapse = ", "),
#            ")", sep = "")
#        if (argname == "...") sel else paste(argname, "=", sel)
#    } else ""
#        res
#    }
#    assign("result", result, envir = resenv)
#    select <- function() tkfocus(tl)
#    assign("select", select, envir = resenv)
#    return(resenv)
#}
#
#guiDlg <- function (title = "Input", message = NULL, help = NULL, sep = NULL,
#width = 50, labelwidth = 0, panes = list(list(type = "entry",
#message = "Enter data:", default = NULL)), GUI = getOption("guiWidgets"))
#{
#    ## Compute a guiDlg object
#    res <- list(list(title = title, message = message, help = help, sep = sep,
#        width = width, labelwidth = labelwidth))
#    ## Add panes
#    if (!is.null(panes) && length(panes) > 0)
#        for (i in 1:length(panes))
#            res[[i + 1]] <- panes[[i]]
#    class(res) <- c("guiDlg", "gui")
#    return(res)
#}
#
#guiDlgFunction <- function (fun, template = NULL, maxargs = 7, var = "res",
#width = 40, labelwidth = 10, displayit = TRUE, execfun = getOption("guiExecFun"))
#{
#    ## This dialog prompts for arguments, given a function
#    ## and it constructs the corresponding command
#    ## fun is the name of a function
#    ## template is an alternate template
#    ## displayit displays the dialog box and get results
#    ## execfun is the function to call to run it
#
#    ## Get fun
#    if (!exists(fun, where = 1, mode = "function"))
#        stop(fun, "does not exist or is not a function!")
#    ## Get formal arguments for this function
#    Form <- formals(get(fun, pos = 1, mode = "function"))
#### TODO: use an existing template
#### TODO: deal with S3 and S4 generic functions!
#### TODO: use syntax for call arg by position!
#    ## Construct a default template for this function
#    if (isHelp(fun)["help"]) {
#        hlp <- function (...) help(...)  # To avoid warning on R CMD check!
#        ## help() function is changed in R 2.10!
#        if (exists("getRversion", mode = "function") &&
#            getRversion() >= '2.10') {
#            Help <- paste("browseURL('", hlp(fun, help_type = "html"), "')",
#                sep = "")
#        } else {  # This is R <= 2.9.x
#            Help <- paste("browseURL('", hlp(fun, htmlhelp = TRUE), "')",
#                sep = "")
#            ## Or simply use: paste("help('", fun, "')", sep = "")
#            ## to use default help system
#        }
#    } else Help <- NULL
#    Tpl <- list(list(fun = fun, var = var, title = "Function assistant",
#        message = NULL, help = Help, sep = NULL, width = width,
#        labelwidth = labelwidth))
#    if (!is.null(Form)) {  # If there are arguments
#        Nargs <- length(Form)
#        ArgsNames <- names(Form)
#        ## Take at most maxargs argument (if more, the rest is ...)
#        if (Nargs > maxargs) N <- maxargs else N <- Nargs
#        ## Create an entry in the list for each arg
#        for (i in 1:N)
#            Tpl[[i + 1]] <- list(type = "entry", argname = ArgsNames[i],
#        default = deparse(Form[[i]]))
#        if (Nargs > maxargs) {  # Include "..."
#            ## Process a message with other args
#            ArgsNames <- ArgsNames[-(1:N)]
#            Form <- Form[-(1:N)]
#### TODO: use deparse here also!
#            OtherArgs <- paste("Other arguments:", paste(ArgsNames, Form,
#                sep = " = ", collapse = ", "))
#            Tpl[[maxargs + 2]] <- list(argname = "...", message = OtherArgs,
#                fixedfont = TRUE)
#        }
#    }
#    class(Tpl) <- c("guiDlg", "gui")
#
#    ## Do we have to return this template or to run it?
#    if (!displayit) return(Tpl)
#
#    ## Otherwise display the dialog box... and get results
#    res <- display(Tpl)
#    ## Do we have to execute it?
#    if (is.null(execfun)) execfun <- "guiEval" # Default evaluator
#    if (execfun != "") {
#		if (exists(execfun, where = -1, mode = "function")) {
#        	get(execfun, pos = -1, mode = "function")(res)
#		} else warning(execfun, " not found!")
# 	}
#    ## Return res invisibly
#    return(invisible(res))
#}
#
#guiEval <- function (code, ident = "GUI ")
#{
#    ## This function is used by default to evaluate constructed code
#    if (is.null(code) || is.na(code) || !inherits(code, "character") ||
#        length(code) == 0) return()
#    ## Echo command
#    Prompt <- getOption("prompt")
#    if (ident != "")
#        Prompt <- paste(Prompt, ident, Prompt, sep = "")
#    cat(Prompt, code[1], "\n", sep = "")
#    if (length(code) > 1) {
#        Continue <- getOption("continue")
#        for ( i in 2:length(code))
#            cat(Continue, code[i], "\n", sep = "")
#    }
#    ## Evaluate this command
#    e <- try(parse(text = code))
#    if (inherits(e, "try-error"))
#        stop("Syntax error!")
#    yy <- withVisible(eval(e, envir = .GlobalEnv))
#    if (yy$visible) print(yy$value)
#}
#
#display <- function (x, ...)
#    UseMethod("display")
#
#display.guiDlg <- function (x, parent = 0, GUI = getOption("guiWidgets"),
#debug = FALSE, ...)
#{
#    ## Check arguments
#    if (!inherits(x, "guiDlg"))
#        stop("'x' must be a guiDlg object!")
#### TODO: check parent
#    if (!is.null(debug) && !is.na(debug)) debug <- (debug == TRUE) else
#        debug <- FALSE
#    if (!inherits(GUI, "character") && !is.null(GUI))
#        stop("'GUI' must be a character string or NULL!")
#
#    ## Do we need to use a different widget than Tcl/Tk?
#    if (!is.null(GUI) && GUI != "tcltk") {  # Custom GUI widgets
#        ## Look for a display.guiDlg.<GUI> function
#        fun <- paste("display.guiDlg", GUI, sep=".")
#        if (exists(fun, where = 1, mode = "function", inherits = TRUE)) {
#            res <- get(fun, pos = 1, mode = "function", inherits = TRUE)(
#            x = x, parent = parent, debug = debug)
#            if (!is.null(res)) {
#                return(res)
#            } else warning("Using default Tcl/tk dialog box instead!")
#        }
#    }
#
#    ## Otherwise, use the default Tcl/Tk dialog box
#    ## Check the content of 'x'
#    X <- x[[1]]
#    panes <- x
#    panes[[1]] <- NULL
#    if (!inherits(X$title, "character") && length(X$title) < 1)
#        stop("'title' must be a non empty character string!")
#    title <- X$title[1]  # Keep only first item for title
#    if (!is.null(X$message)) message <- paste(as.character(X$message),
#        collapse = "\n") else message <- NULL
#    if (!is.null(X$help) && !inherits(X$help, "character"))
#        stop("'help' must be NULL or a character string!")
#    help <- X$help[1]  # Keep only first item
#	if (is.null(X$labelwidth)) X$labelwidth <- 0  # Default value
#    if (!is.numeric(X$labelwidth))
#        stop("'labelwidth' must be a number or NULL!")
#    labelwidth <- as.integer(X$labelwidth)
#    if (labelwidth < 0) labelwidth <- 0
#    if (is.null(X$width)) X$width <- 40  # Default value
#    if (!is.numeric(X$width))
#        stop("'width' must be a number or NULL!")
#    width <- as.integer(X$width)
#    ## If "Help" is displayed min(width + labelwidth) = 35 else it is 20
#    if (is.null(help)) minwidth <- 20 - labelwidth else minwidth <- 35 - labelwidth
#    if (width < minwidth) width <- minwidth
#    if (width < 10) width <- 10 # Minimum absolute width of 10
#### TODO: check these arguments
#    fun <- X$fun
#    var <- X$var
#
#    ## Make sure style is defined
#    style <- guiSetStyle.tcltk(getOption("guiStyle"))
#    ## Size widgets according to text font measure
#    pixwidth <- as.integer((width + labelwidth) * style$font.measure["text"])
#    ## Do we need to use pane separators?
#    if (is.null(X$sep)) sep <- style$sep else sep <- (X$sep == TRUE)
#
#    ## Construct the dialog box
#    cnt <- tktoplevel(class = "guiDlg")
#    tkwm.withdraw(cnt)  # Do not show it until it is completelly constructed!
#    on.exit(tkdestroy(cnt))  # Make sure we don't left it open in case of error!
#    tktitle(cnt) <- title
#    ## Do we need to display a "header" for a function construction?
#    banner <- FALSE
#    txtAssign <- NULL
#    if (!is.null(fun)) {
#        banner <- TRUE
#        funFrame <- tkframe(cnt)
#        if (!is.null(var)) {  # Give the possibility to assign to a variable
#            varAssign <- tclVar(var)
#            txtAssign <- tkentry(funFrame, textvariable = varAssign,
#                width = max(labelwidth, 10), font = style$font.text, fg = style$fg[1],
#                background = "white", relief = style$relief)
#            #tkconfigure(txtAssign, validate = "key", validatecommand = onVarEdit)
#            tkselection.from(txtAssign, "0")
#            tkselection.to(txtAssign, "end")
#            tkicursor(txtAssign, "end")
#            labFun <- tklabel(funFrame, text = paste("<- ", fun, "()", sep = ""),
#                font = style$font.fixed, justify = "left", fg = style$fg[4])
#            tkgrid(txtAssign, labFun)
#            tkgrid(funFrame, sticky = "w", padx = style$pads[1], pady = style$pads[3])
#        } else {  # No assignation allowed
#            labFun <- tklabel(funFrame, text = paste(fun, "()", sep = ""),
#                font = style$font.fixed, justify = "left", fg = style$fg[4])
#            if (labelwidth > 0) {
#                labSpacer <- tklabel(funFrame, text = " ", font = style$font.fixed,
#                    width = labelwidth)
#                tkgrid(labSpacer, labFun)
#            } else tkgrid(labFun)
#            tkgrid(funFrame, sticky = "w", padx = style$pads[1], pady = style$pads[3])
#        }
#    }
#    if (!is.null(message)) {  # Display a banner with the message
#        banner <- TRUE
#        dialoglabel <- tklabel(cnt, text = message, font = style$font.emph,
#            justify = "left", fg = style$fg[5], wraplength = pixwidth)
#        tkgrid(dialoglabel, sticky = "w", padx = style$pads[1])
#    }
#    if (sep && banner) {
#        sepa <- tkcanvas(cnt, height = "0",relief = "groove", borderwidth = "1",
#            width = pixwidth)
#        tkgrid(sepa)
#    }
#
#    ## Construct a dlg object
#    dlg <- list(call = match.call(), title = title, message = message,
#        container = cnt, style = style, width = width, pixwidth = pixwidth,
#        labelwidth = labelwidth, sep = sep, panes = panes, result = character(0))
#    class(dlg) <- c("guiDlg", "gui")
#
#    ## Call guiPane.tcltk() to construct the panes
#    for (i in 1:length(panes)) {
#        dlg <- guiPane.tcltk(dlg, i)
#        if (sep) {
#            sepa <- tkcanvas(cnt, height = "0",relief = "groove", borderwidth = "1",
#                width = pixwidth)
#            tkgrid(sepa)
#        }
#    }
#
#    ## Since onOk must update dlg$result, but I cannot pass dlg
#    ## as argument to and from the onOk function, I save it in a temporary variable
#    vardlg <- tempvar(".dlg")  # Needed to store state of the dialog box
#    assign(vardlg, dlg, pos = 1)
#    if (!debug) on.exit(remove(list = vardlg, pos = 1), add = TRUE)
#    getdlg <- eval(parse(text = paste("function() get('", vardlg, "', pos = 1)",
#        sep = "")))
#    setdlg <- eval(parse(text = paste("function(dlg) ", vardlg, " <<- dlg", sep = "")))
#    onOk <- function () {
#        dlg <- getdlg()  # Retrieve the dialog object from the temp variable
#        panes <- dlg$panes
#        ## Get results from individual panes
#        res <- NULL
#        for (i in 1:length(panes))
#            res[i] <- eval(parse(text = "result()"), envir = panes[[i]]$env)
#        dlg$result <- res
#        setdlg(dlg)
#		## Indicate we clicked 'OK'
#		assignTemp(".guiDialog.res", "ok")
#        tkdestroy(cnt)
#    }
#    onCancel <- function() tkdestroy(cnt)
#
#    ## Add the dialog buttons
#    butFrame <- tkframe(cnt)
#    butOk <- tkbutton(butFrame, text = "OK", width = "10", command = onOk,
#        default = "active", font = style$font.label, fg = style$fg[2])
#    labSep <- tklabel(butFrame, text = " ", font = style$font.label)
#    butCancel <- tkbutton(butFrame, text = "Cancel", width = "10",
#        command = onCancel, font = style$font.label, fg = style$fg[3])
#    if (is.null(help)) {
#        tkgrid(butOk, labSep, butCancel, sticky = "w")
#    } else {
#        labSep2 <- tklabel(butFrame, text = "     ", font = style$font.label)
#        onHelp <- function() eval(parse(text = help), envir = .GlobalEnv)
#        butHelp <- tkbutton(butFrame, text = "Help", width ="10",
#            command = onHelp, font = style$font.label, fg = style$fg[1])
#        tkgrid(butOk, labSep, butCancel, labSep2, butHelp, sticky = "w")
#        tkbind(cnt, "<F1>", onHelp)
#    }
#    tkgrid(butFrame, padx = style$pads[1], pady = style$pads[2])
#    ## Finalize the configuration of the dialog box
#    tkwm.resizable(cnt, 0, 0)
#    tkwm.protocol(cnt, "WM_DELETE_WINDOW", onCancel)
#    tkbind(cnt, "<Return>", onOk)
#    tkbind(cnt, "<Escape>", onCancel)
#    ## The only solution I have found to eliminate minbutton and make the dialog
#    ## box always on top of R Console under Windows is the following one (to change!)
#    if (.Platform$OS.type == "windows")
#        tcl("wm", "attributes", cnt, toolwindow = 1, topmost = 1)
#    .Tcl("update idletasks")
#    tkwm.deiconify(cnt)
#    ## tkwm.deiconify() is enough! tkfocus(force = cnt)
#    tkgrab.set(cnt)# This is a modal dialog box => keep focus!
#    if (is.null(txtAssign)) {
#        ## Select adequate widget in first pane
#        eval(parse(text = "select()"), envir = dlg$panes[[1]]$env)
#    } else tkfocus(txtAssign)
#    ## Set by default return value of the dialog box to "cancel"
#    assignTemp(".guiDialog.res", "cancel")
#    tkwait.window(cnt)
#    ## Did we cancelled the dialog box?
#	if (get(".guiDialog.res", envir = .GlobalEnv) == "cancel") return(NULL)
#	## Get the updated version of the dialog box
#    dlg <- get(vardlg, pos = 1)
#    res <- dlg$result
#    ## If this is a function, compute corresponding R code
#    if (!is.null(fun)) {
#        res <- paste(res[res != ""], collapse = ", ")
#        res <- paste(fun, "(", res, ")", sep = "")
#        if (!is.null(txtAssign) && (varname <- tclvalue(varAssign)) != "")
#            res <- (paste(varname, "<-", res))
#        res <- strwrap(res, exdent = 4)
#    }
#    return(res)
#}
