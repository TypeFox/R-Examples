.koUnit <- function (cmd, warn = FALSE, ...)
{
	## Look if koCmd() exists, otherwise, we are probably not connected to Komodo
	if (exists("koCmd", mode = "function")) {
		res <- get("koCmd")(cmd, ...)
		if (isTRUE(warn) & inherits(res, "try-error"))
			warning("Komodo is not available or did not process this command correctly")
		return(res)
	} else {
		if (isTRUE(warn))
			warning("You must establish a connection with Komodo/SciViews-K to use this function")
	}
}

koUnit_isAutoTest <- function ()
{
	res <- .koUnit('sv.socket.serverWrite(sv.r.unit.isAutoTest());')
	return(res == "true")
}

koUnit_setAutoTest <- function (state)
{
	if (isTRUE(state)) state <- "true" else state <- "false"
	res <- .koUnit('sv.r.unit.setAutoTest(<<<data>>>);', data = state)
}

koUnit_runTest <- function ()
	res <- .koUnit('sv.r.unit.runTest();')

koUnit_showRUnitPane <- function (state)
{
	if (missing(state)) state <- ""
	else if (isTRUE(state)) state <- "true" else state <- "false"
	res <- .koUnit('sv.r.unit.showRUnitPane(<<<data>>>);', data = state)
}

koUnit_version <- function ()
{
	res <- .koUnit('sv.socket.serverWrite(sv.r.unit.version + "." + sv.r.unit.release);')
	return(res)
}
