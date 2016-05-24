print.guiWin <- function (x, ...)
{
	cat("A SciViews GUI window object:", "\n")
	print(unclass(x))
	invisible(x)
}

winAdd <- function (name = "win1", type = "tkWin", parent = .TkRoot,
title = NULL, pos = NULL, bind.delete = TRUE, ...)
{
	## Add a window. This mechanism should be able to use different kinds of
	## graphical widgets, but currently, only Tcl/Tk is supported.
	invisible(switch(type,
		tkWin = tkWinAdd(name = name, parent = parent, title = title, pos = pos,
			bind.delete = bind.delete, ...),
		stop("Only type = \"tkWin\" is currently supported")))
}

winDel <- function (window)
{
	## Process depends on the kind of window to delete
	## Currently, only Tk windows are supported
	if (inherits(winGet(window), "tkguiWin")) {
		invisible(tkWinDel(window))
	} else stop("Unsupported window type")
}

winGet <- function (window)
{
	## Retrieve a "guiWin" object from .guiWins, given its name
	return(getTemp(".guiWins")[[window]])
}

winNames <- function ()
{
	## List all recorded windows in .guiWins
### TODO: if Rgui, list also console, graph, editors and pagers!
	res <- names(getTemp(".guiWins"))
	if (is.null(res)) res <- character(0)
	res
}
