## A Tcl/Tk version of the directory selection dialog box
dlgDir.tcltkGUI <- function (default = getwd(), title, ..., gui = .GUI)
{
    gui$setUI(widgets = "tcltkGUI")
	## TODO: we don't display multiline title here! => what to do?
    res <- tclvalue(tkchooseDirectory(initialdir = gui$args$default,
		mustexist = FALSE, title = gui$args$title))
    ## tkchooseDirectory returns "" if cancelled
	if (res == "") res <- character(0)
    gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}
