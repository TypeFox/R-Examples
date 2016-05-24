## Define the S3 method
## TODO: define default extension!!!
dlgSave <- function (default, title, filters = dlgFilters["All", ],
..., gui = .GUI)
{
    ## A 'save file' dialog box
    ## title is used as caption of the dialog box
    ## default allows to preselect a file
    ## Always ask for confirmation in case the file already exists
    ## filters is a n x 2 matrix of characters with description and filter
    ## for instance: "R or S files (*.R, *.q)"       "*.R;*.q"
    ## It could be also an even number of character strings that will be
	## reorganized into a n x 2 matrix.
	
	if (missing(default) || !length(default)) default <- character(0)
	if (!gui$startUI("dlgSave", call = match.call(), default = default,
		msg = "Displaying a modal save file dialog box",
		msg.no.ask = "A modal save file dialog box was by-passed"))
		return(invisible(gui))
	
	## Check and rework main arguments and place them in gui$args
	if (missing(default) || is.null(default))
		default <- file.path(path.expand(getwd()), "untitled", sep = "")
	if (!length(default)) default <- NULL
	if (!is.null(default)) {
		default <- as.character(default)[1]
		## Under Windows, it uses \\ as separator, although .Platform$file.sep
		## is now / (tested in R 2.11.1) => replace it
		if (.Platform$OS.type == "windows")
			default <- gsub("\\\\", "/", default)
		## Check that dir of default already exists
		dir <- dirname(default)
		## If not there, or not a dire, replace by current working dir...
		if (!file.exists(dir) || !file.info(dir)$isdir)
			default <- file.path(getwd(), basename(default))
	}
	if (missing(title) || title == "") {
		title <- "Save file as"
	} else title <- as.character(title)[1]
	## Check that filter is a nx2 character matrix, or try reshape it as such
	if (is.matrix(filters)) {
		if (ncol(filters) != 2 || !is.character(filters))
			filters <- NULL
	} else {
		if (length(filters) %% 2 != 0) {
			filters <- NULL
		} else { # Try to reshape it
			filters <- matrix(as.character(filters), ncol = 2, byrow = TRUE)
		}
	}
	gui$setUI(args = list(default = default, title = title, filters = filters))
	
	## ... and dispatch to the method
	UseMethod("dlgSave", gui)
}

## Used to break the chain of NextMethod(), searching for a usable method
## in the current context
dlgSave.gui <- function (default, title, filters = dlgFilters["All", ], ..., gui = .GUI)
{
	msg <- paste("No workable method available to display a file save dialog box using:",
		paste(guiWidgets(gui), collapse = ", "))
	gui$setUI(status = "error", msg = msg, widgets = "none")
	stop(msg)
}

## The pure textual version used as fallback in case no GUI could be used
dlgSave.textCLI <- function (default, title, filters = dlgFilters["All", ],
..., gui = .GUI)
{
	gui$setUI(widgets = "textCLI")
	## Ask for the file
	res <- readline(paste(gui$args$title, " [", gui$args$default, "]: ",
		sep = ""))
	if (res == "") res <- gui$args$default else res <- res
	## In case we pasted a string with single, or double quotes, or spaces
	## eliminate them
	res <- sub("^['\" ]+", "", sub("['\" ]+$", "", res))
	if (length(res)) {
		res <- normalizePath(res)
		## If file already exists => ask for confirmation...
		if (file.exists(res)) {
			choices <- c("ok", "cancel")
			ret <- select.list(choices,
				title = "Confirm you want to replace this file", graphics = FALSE)
			if (ret == "" || ret == "cancel") res <- character(0) # Cancelled
		}
	}
	gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}

## The native version of the file save dialog box
dlgSave.nativeGUI <- function (default, title, filters = dlgFilters["All", ],
..., gui = .GUI)
{
	gui$setUI(widgets = "nativeGUI")
	## A 'save file' dialog box
    ## If cancelled, then return character(0)
    ## This dialog box is always modal
	##
	## It is a replacement for choose.files(), tkgetSaveFile()
	## & file.choose(new = TRUE), not implemented yet in R 2.14, by the way
	res <- switch(Sys.info()["sysname"],
		Windows = .winDlgSave(gui$args$default, gui$args$title, gui$args$filters),
		Darwin = .macDlgSave(gui$args$default, gui$args$title, gui$args$filters),
		.unixDlgSave(gui$args$default, gui$args$title, gui$args$filters)
	)
	
	## Do we need to further dispatch?
	if (is.null(res)) NextMethod("dlgSave", gui) else {
		gui$setUI(res = res, status = NULL)
		return(invisible(gui))
	}
}

## Windows version
.winDlgSave <- function (default, title, filters = dlgFilters["All", ])
{
	## Note: this dialog box is a very bad design for saving as dialog box:
	## It displays the "Open" button and complains if the file does not exist!
	## TODO: should be replaced by the Tk version for now
	if (!is.matrix(filters)) filters <- matrix(filters, ncol = 2, byrow = TRUE)
	if (!length(default)) default <- ""
	res <- choose.files(default = default, caption = title,
		multi = FALSE, filters = filters, index = 1)
    if (length(res)) {
		res <-  gsub("\\\\", "/", res)
		if (file.exists(res) && tolower(winDialog(type = "okcancel",
			message = "The file already exists. It will be replaced!")) == "cancel")
			res <- character(0) # User cancelled!
	}
	return(res)
}

## Mac OS X version
.macDlgSave <- function (default, title, filters = dlgFilters["All", ])
{
    ## TODO: filters are implemented differently on the Mac => how to do this???
	if (!is.matrix(filters)) filters <- matrix(filters, ncol = 2, byrow = TRUE)
	## Display a modal file save selector with native Mac dialog box
	#if (.Platform$GUI == "AQUA") app <- "(name of application \"R\")" else
	## This works from Mac OS X 10.5 Leopard:
	if (.Platform$GUI == "AQUA") app <- "(name of application id \"Rgui\")" else
		if (.isJGR()) app <- "\"JGR\"" else app <- "\"Terminal\""
	## Avoid displaying warning message when the user clicks on 'Cancel'
	owarn <- getOption("warn")
	on.exit(options(warn = owarn))
	options(warn = -1)
	if (title == "") mcmd <- "" else mcmd <- paste("with prompt \"",
		title, "\" ", sep = "")
	if (length(default) && default != "") {
		## Default dir must be an existing dir... otherwise, the cmd fails!
		defdir <- dirname(default)
		if (!file.exists(defdir) || !file.info(defdir)$isdir)
			defdir <- getwd()
		mcmd <- paste(mcmd, " default location \"", defdir, "\"", sep = "")
		deffile <- basename(default)
		if (deffile != "*.*" && deffile != "*")
			mcmd <- paste(mcmd, " default name \"", deffile, "\"", sep = "")
	}
	cmd <- paste("-e 'tell application ", app,
		" to set filename to choose file name ", mcmd,
		"' -e 'POSIX path of filename'", sep = "")
	## I cannot use system(intern = TRUE) with this in R.app/R64.app
	## (deadlock situation?), but I can in R run in a terminal. system2() also
	## works, but this preclue of using svDialogs on R < 2.12.0.
	## The hack is thus to redirect output to a file, then, to read the content
	## of that file and to destroy it
	tfile <- tempfile()
	on.exit(unlink(tfile))
	res <- try(system(paste("osascript", cmd, ">", tfile), wait = TRUE,
		intern = FALSE, ignore.stderr = TRUE), silent = TRUE)
	if (inherits(res, "try-error") || !length(res)) return(character(0))
	if (res > 0) return(character(0)) # User cancelled input
	res <- readLines(tfile)
	res <- res[res != ""] # Eliminate empty lines
	## Note: confirmation of replacement is built-in here
	return(res)	
}

## Linux/Unix version
## TODO: if no extension provided, displays '.' => make sure to change this!
.unixDlgSave <- function (default, title, filters = dlgFilters["All", ])
{
    ## Note: only existing filenames can be selected as default, otherwise, the
	## argument is ignored!
	## zenity must be installed on this machine!
    if (Sys.which("zenity") == "") return(NULL)
    ## Avoid displaying warning message in case user clicks on Cancel
    owarn <- getOption("warn")
    on.exit(options(warn = owarn))
    options(warn = -1)
    ## Use zenity to display the file save selection
    ## Construct the -file-filter options
	fcmd <- ""
	if (!is.matrix(filters)) filters <- matrix(filters, ncol = 2, byrow = TRUE)
	nf <- nrow(filters)
	if (nf > 0) for (i in 1:nf)
		fcmd <- paste(fcmd, " --file-filter=\"", filters[i, 1], " | ",
			gsub(";", " ", filters[i, 2]), "\"", sep = "")
	msg <- paste("zenity --file-selection --save --title=\"", title,
	"\" --filename=\"", default, "\" ", fcmd, sep = "")
	res <- system(msg, intern = TRUE, ignore.stderr = TRUE) # Because error if file not found!
	if (!length(res)) return(character(0)) else {
		if (file.exists(res)) { # Ask for confirmation!
			msg <- paste("zenity --question --text=\"",
				"This file already exists. It will be replaced!",
				"\" --ok-label=\"OK\" --cancel-label=\"Cancel\"",
				" --title=\"Question\"", sep = "") 
			if (system(msg) > 0) return(character(0)) # Cancelled, or another error
		}
		return(strsplit(res, "|", fixed = TRUE)[[1]])
	}
}
