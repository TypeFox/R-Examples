## A Tcl/Tk version of the file save dialog box
dlgSave.tcltkGUI <- function (default, title, filters = dlgFilters["All", ], ..., gui = .GUI)
{
    gui$setUI(widgets = "tcltkGUI")
    ## In tkgetSaveFile, filters are presented differently!
    filters <- gui$args$filters
	filters <- paste("{\"", filters[, 1], "\" {\"", gsub(";", "\" \"",
		filters[, 2]), "\"}}", sep = "", collapse = " ")
    ## Use tkgetSaveFile()
	default <- gui$args$default
	if (!length(default) || default == "")
		default <- file.path(getwd(), "*.*")
	res <- as.character(tkgetSaveFile(title = gui$args$title,
		initialfile = basename(default),
        initialdir = dirname(default), # defaultextension = Not used!
		filetypes = filters))
	if (!length(res) || res == "") res <- character(0)
	## Note: confirmation is built-in!
    gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}
