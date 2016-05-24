setUI <- function (..., gui = .GUI)
	UseMethod("setUI", gui)

setUI.gui <- function (fun, call, args, res, widgets, status, msg = NULL,
..., gui = .GUI)
{
	## Set something in the GUI
	if (!missing(call)) {
		if (!inherits(call, "call"))
			stop("'call' must be a call expression (use match.call)")
		if (missing(fun))
			stop("'fun' must be provided with call (original name of the calling function)")
		fun <- as.character(fun)[1]
		## Rework call to make sure to have original name for function and gui
		call[1] <- call(fun)
		if (!is.null(gui$name)) call["gui"] <- call(as.character(gui$name)[1])
		gui$call <- call
	}
	if (!missing(args)) gui$args <- args
	if (!missing(res)) gui$res <- res
	if (!missing(widgets)) gui$widgets <- widgets
	if (!missing(status)) {
		comment(status) <- msg
		gui$status <- status
	} else if (!is.null(msg))
		stop("You must provide 'status' at the same time as 'msg'")
	## Process ...
	moreArgs <- list(...)
	l <- length(moreArgs)
	n <- names(moreArgs)
	if (l) for (i in 1:l) {
		gui[[n[i]]] <- moreArgs[i]
	}
	return(invisible(gui))
}

startUI <- function (..., gui = .GUI)
	UseMethod("startUI", gui)

startUI.gui <- function (fun, call, default, widgets = NULL,
status = "busy-modal", msg = "Displaying a modal dialog box",
msg.no.ask = "A modal dialog box was by-passed", ..., gui = .GUI)
{
	if (dontAsk(gui)) { ## Don't ask, just return default value
		gui$setUI(fun = fun, call = call, res = default, 
			widgets = "none", status = "by-passed", msg = msg.no.ask)
		return(FALSE)
	} else {
		## Inform that we are about to display a modal dialog box
		gui$setUI(fun = fun, call = call, res = NA,
			widgets = NULL, status = status, msg = msg)
		return(TRUE)
	}
}
