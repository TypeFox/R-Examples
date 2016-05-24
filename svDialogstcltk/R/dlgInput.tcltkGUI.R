## A Tcl/Tk implementation of the input box
dlgInput.tcltkGUI <- function (message = "Enter a value", default = "", ...,
gui = .GUI)
{
    gui$setUI(widgets = "tcltkGUI")
	## A simple text input box
    ## This dialog box is always modal
    ## Return either a string, or character(0) if 'Cancel' clicked
	## Do we have access to ttk widgets?
	if (as.character(tcl("info", "tclversion")) < "8.5") {
        ttklabel <- tklabel
		ttkbutton <- tkbutton
		ttkentry <- tkentry
		ttkframe <- tkframe
	}
	inbox <- tktoplevel()
	tkwm.withdraw(inbox)
    ## Construct the dialog box
	tktitle(inbox) <- "Question"
    ## The message
	wlabel <- ttklabel(inbox, text = gui$args$message)
    tkgrid(wlabel, sticky = "w", padx = 5, pady = 5)
	## The text entry
	varText <- tclVar(gui$args$default)
    wtext <- ttkentry(inbox, width = "50", textvariable = varText)
    tkgrid(wtext, padx = 5, pady = 5)
	## The buttons
	onOk <- function () {
        tkdestroy(inbox)
    }
	onCancel <- function () {
		tkdestroy(inbox)
		## Change the text to @@@@@Cancelled!@@@@@ in varText
		tclvalue(varText) <- "@@@@@Cancelled!@@@@@"
	}
    butFrame <- ttkframe(inbox)
    wbutOK <- ttkbutton(butFrame, text = "OK", width = "12", command = onOk,
        default = "active")
    wbutCancel <- ttkbutton(butFrame, text = "Cancel", width = "12",
        command = onCancel)
    if (.Platform$OS.type == "windows") {
		tkgrid(wbutOK, wbutCancel, sticky = "w", padx = 5)
	} else {
		tkgrid(wbutCancel, wbutOK, sticky = "w", padx = 5)
	}
    tkgrid(butFrame, pady = 5)
    ## Other event bindings
	tkwm.protocol(inbox, "WM_DELETE_WINDOW", onCancel)
    tkbind(inbox, "<Return>", onOk)
	tkbind(inbox, "<KeyPress-Escape>", onCancel)
    .Tcl("update idletasks")
    tkselection.range(wtext, "0", "100000")
    tkicursor(wtext, as.character(nchar(gui$args$default)))
	## Center the window on screen
	W <- as.integer(tclvalue(tcl("winfo", "screenwidth", inbox)))
	H <- as.integer(tclvalue(tcl("winfo", "screenheight", inbox)))
	X <- as.integer(tclvalue(tcl("winfo", "width", inbox)))
	Y <- as.integer(tclvalue(tcl("winfo", "height", inbox)))
	X0 <- W/2 - X/2 - as.integer(tclvalue(tcl("winfo", "vrootx", inbox)))
	Y0 <- H/2 - Y/2 - as.integer(tclvalue(tcl("winfo", "vrooty", inbox)))
	pos <- paste("+", round(X0), "+", round(Y0), sep = "")
	tkwm.geometry(inbox, pos)
	## Do noy allow to resize the window
	tkwm.resizable(inbox, 0, 0)
	## Show it now and wait for it is destroyed
	tkwm.deiconify(inbox)
	## Make sure the window is always on top
	try(tcl("wm", "attributes", inbox, topmost = 1), silent = TRUE)
	## Focus now inside the entry widget
	tkfocus(wtext)
	## Wait for the window to be destroyed
	tkwait.window(inbox)
	## Get the value
	res <- tclvalue(varText)
	## Check if the user cancelled the dialog box
	if (res == "@@@@@Cancelled!@@@@@") res <- character(0)
    gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}
