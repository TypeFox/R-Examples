## Because guiChange() can manage all cases where the GUI already exists or not
## and guiAdd() must also manage the case where the GUI already exists,
## the two functions are pretty much exchangeable!
guiAdd <- function (gui.name = ".GUI", widgets = c("nativeGUI", "textCLI"), ask)
	guiChange(gui.name = gui.name, widgets = widgets, ask = ask, reset = FALSE)

guiChange <- function (gui.name = ".GUI", widgets = c("nativeGUI", "textCLI"),
reset = FALSE, ask)
{
	gui.name <- as.character(gui.name)[1]
	## Do the object already exists in SciViews:TempEnv?
	if (exists(gui.name, envir = .TempEnv(), inherits = FALSE)) {
		guiObj <- get(gui.name, envir = .TempEnv(), inherits = FALSE)
		if (!is.gui(guiObj))
			stop("'gui.name' must be a character string naming a 'gui' object in SciViews:TempEnv")
		guiWidgets(guiObj, reset = reset) <- widgets
		if (isTRUE(reset)) {
			## Make sure name is correct
			guiObj$name <- gui.name
			## Use default for ask, if not provided
			if (missing(ask)) guiObj$ask <- NULL
		}
	} else { # Create a new 'gui' object in SciViews:TempEnv
		if (is.na(gui.name))
			stop("Wrong 'gui.name', provide a length 1 character string")
		guiObj <- new.env(parent = .GlobalEnv)
		guiObj$name <- gui.name
		guiObj$ask <- NULL
		guiObj$call <- NULL
		guiObj$res <- NULL
		guiObj$status <- NULL
		guiObj$widgets <- NULL
		## Define initial inheritance
		class(guiObj) <- unique(c(widgets, "gui", "environment"))
		assign(gui.name, guiObj, envir = .TempEnv())
	}
	## Do we change the ask value?
	if (!missing(ask))
		if (is.null(ask)) guiObj$ask <- NULL else
			guiObj$ask <- isTRUE(as.logical(ask))
	## Return the modified GUI
	return(guiObj)
}

guiRemove <- function (gui.name)
{
	## Eliminate the corresponding variable, after some housekeeping
	if (gui.name == ".GUI")
		stop("You cannot delete the default GUI named '.GUI'! Maybe use ?guiChange.")
	if (!exists(gui.name, envir = .TempEnv(), inherits = FALSE))
		return(invisible(FALSE))
	## TODO: housekeeping here!
	rm(list = gui.name, envir = .TempEnv(), inherits = FALSE)
	return(invisible(TRUE))	
}

## List all GUI objects found in SciViews:TempEnv
guiList <- function ()
{
	lst <- ls(envir = .TempEnv(), all.names = TRUE)
	if (!length(lst)) return(character(0))
	## Check which item inherits from 'gui'
	res <- lst[sapply(lst, function (x) is.gui(get(x, envir = .TempEnv(),
		inherits = FALSE)))]
	return(res)
}

guiWidgets <- function (gui, gui.name = ".GUI")
{
	if (missing(gui))
		if (exists(gui.name, envir = .TempEnv(), inherits = FALSE)) {
			gui <- get(gui.name, envir = .TempEnv(), inherits = FALSE)
		} else stop("'gui' object '", gui.name, "' not found")
	if (!is.gui(gui)) stop("Provide a 'gui' object or its name")
	classes <- class(gui)
	## Keep only all classes before 'gui'
	classes <- classes[!cumsum(classes == "gui")]
	return(classes)
}

`guiWidgets<-` <- function (x, reset = FALSE, value)
{
	value <- as.character(value)
	if (isTRUE(as.logical(reset))) {
		## Change completely class
		class(x) <- unique(c(value, "gui", "environment"))
	} else {
		## Add 'value' items to current classes
		classes <- class(x)
		value <- value[!value %in% classes]
		if (length(value))
			class(x) <- c(value, classes)
	}
	return(x)
}

guiAsk <- function (gui.or.name, ask)
{
	if (missing(gui.or.name)) {
		## Query or change the default value in 'gui.ask' option
		if (missing(ask)) {
			res <- getOption("gui.ask", default = NULL)
			if (is.null(res)) res <- structure(TRUE, comment = "default")
			return(res)
		}
		if (!is.null(ask)) ask <- isTRUE(as.logical(ask))
		res <- options(gui.ask = ask)$gui.ask
		if (!is.null(res)) res <- as.logical(res)
	} else {
		if (is.gui(gui.or.name)) guiObj <- gui.or.name else {
			if (!exists(gui.or.name, envir = .TempEnv(), inherits = FALSE))
				stop("'gui' object '", gui.or.name, "' not found")
			guiObj <- get(gui.or.name, envir = .TempEnv(), inherits = FALSE)
			if (!is.gui(guiObj))
				stop("'gui.or.name' must be a 'gui' object in SciViews:TempEnv or its name")
		}

		## Query or set guiObj$ask
		if (missing(ask)) {
			res <- guiObj$ask
			if (is.null(res)) ## Look at default value
				res <- structure(getOption("gui.ask", default = TRUE),
					comment = "default")
			return(res)
		}
		## Change the value for this GUI
		res <- guiObj$ask
		if (is.null(ask)) guiObj$ask <- NULL else
			guiObj$ask <- isTRUE(as.logical(ask))
	}
	## Return default value
	return(invisible(res))
}

`guiAsk<-` <- function (x, value)
{
	if (!is.gui(x))
		stop("guiAsk must be applied to a 'gui' object")
	if (is.null(value)) x$ask <- NULL else
		x$ask <- isTRUE(as.logical(value))
	return(x)
}
