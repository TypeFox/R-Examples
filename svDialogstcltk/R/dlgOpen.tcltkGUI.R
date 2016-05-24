## A Tcl/Tk version of the file open dialog box
dlgOpen.tcltkGUI <- function (default, title, multiple = FALSE,
filters = dlgFilters["All", ], ..., gui = .GUI)
{
    gui$setUI(widgets = "tcltkGUI")
    ## In tkgetOpenFile, filters are presented differently!
    filters <- gui$args$filters
	filters <- paste("{\"", filters[, 1], "\" {\"", gsub(";", "\" \"",
		filters[, 2]), "\"}}", sep = "", collapse = " ")
    ## Use tkgetOpenFile()
    res <- as.character(tkgetOpenFile(title = gui$args$title,
		initialfile = basename(gui$args$default),
        initialdir = dirname(gui$args$default), multiple = gui$args$multiple,
		filetypes = filters))
	if (length(res) == 1 && res == "") res <- character(0)
    gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}
