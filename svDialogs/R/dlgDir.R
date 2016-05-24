## Define the S3 method
dlgDir <- function (default = getwd(), title, ..., gui = .GUI)
{
	if (!gui$startUI("dlgDir", call = match.call(), default = default,
		msg = "Displaying a modal dir selection dialog box",
		msg.no.ask = "A modal dir selection dialog box was by-passed"))
		return(invisible(gui))
	
	## Check and rework main arguments and place them in gui$args
	if (!length(default)) default <- getwd() else
		default <- as.character(default)[1]
    if (file.exists(default))
        if (!file.info(default)$isdir) default <- dirname(default) # Need a dir 
	default <- path.expand(default)
	## Under Windows, it uses \\ as separator, although .Platform$file.sep
	## is now / (tested in R 2.11.1) => replace it
	if (.Platform$OS.type == "windows")
		default <- gsub("\\\\", "/", default)
	if (missing(title) || !length(title) || title == "") {
		title <- "Choose a directory"
	} else {
		title <- paste(as.character(title), collapse = "\n")
	}
	gui$setUI(args = list(default = default, title = title))
	
	## ... and dispatch to the method
	UseMethod("dlgDir", gui)
}

## Used to break the chain of NextMethod(), searching for a usable method
## in the current context
dlgDir.gui <- function (default = getwd(), title, ..., gui = .GUI)
{
	msg <- paste("No workable method available to display a dir selection dialog box using:",
		paste(guiWidgets(gui), collapse = ", "))
	gui$setUI(status = "error", msg = msg, widgets = "none")
	stop(msg)
}

## The pure textual version used a fallback in case no GUI could be used
dlgDir.textCLI <- function (default = getwd(), title, ..., gui = .GUI)
{
	gui$setUI(widgets = "textCLI")
	## Ask for the directory
	res <- readline(paste(gui$args$title, " [", gui$args$default, "]: ",
		sep = ""))
	if (res == "") res <- gui$args$default else res <- res
	## In case we pasted a string with single, or double quotes, or spaces
	## eliminate them
	res <- sub("^['\" ]+", "", sub("['\" ]+$", "", res))
	## To get the same behaviour as the GUI equivalents, we must make sure
	## it is a directory, or try to create it (possibly recursively, if it
	## does not exist). Also return absolute path
	if (file.exists(res)) {
		## Check that this is a directory, not a file!
		if (!file.info(res)$isdir) {
			warning(res, " is not a directory")
			res <- character(0) # Same as if the user did cancel the dialog box
		}
	} else {
		## The directory does not exists, try to create it now...
		dir.create(res, recursive = TRUE)
		if (!file.exists(res) || !file.info(res)$isdir) {
			warning("Error while creating the directory ", res)
			res <- character(0)
		}
	}
	if (length(res)) res <- normalizePath(res)
	gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}

## The native version of the input box
dlgDir.nativeGUI <- function (default = getwd(), title, ..., gui = .GUI)
{
	gui$setUI(widgets = "nativeGUI")
	## A 'choose a directory' dialog box
    ## It almost follows the conventions of tkchooseDirectory()
    ## The argument default indicates the initial directory
    ## If cancelled, then return character(0)
    ## This dialog box is always modal
	##
	## Replacement for choose.dir(), tk_choose.dir() & tkchooseDirectory()
	res <- switch(Sys.info()["sysname"],
		Windows = .winDlgDir(gui$args$default, gui$args$title),
		Darwin = .macDlgDir(gui$args$default, gui$args$title),
		.unixDlgDir(gui$args$default, gui$args$title)
	)
	
	## Do we need to further dispatch?
	if (is.null(res)) NextMethod("dlgDir", gui) else {
		gui$setUI(res = res, status = NULL)
		return(invisible(gui))
	}
}

## Windows version
.winDlgDir <- function (default = getwd(), title = "")
{
	res <- choose.dir(default = default, caption = title)
    if (is.na(res)) res <- character(0) else res <-  gsub("\\\\", "/", res)
	return(res)
}

## Mac OS X version
.macDlgDir <- function (default = getwd(), title = "")
{
    ## Display a modal directory selector with native Mac dialog box
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
	cmd <- paste("-e 'tell application ", app,
		" to set foldername to choose folder ", mcmd, "default location \"",
		default , "\"' -e 'POSIX path of foldername'", sep = "")
	## I cannot use system(intern = TRUE) with this in R.app/R64.app
	## (deadlock situation?), but I can in R run in a terminal. system2() also
	## works, but this preclue of using svDialogs on R < 2.12.0.
	## The hack is thus to redirect output to a file, then, to read the content
	## of that file and to desctroy it
	tfile <- tempfile()
	on.exit(unlink(tfile))
	res <- try(system(paste("osascript", cmd, ">", tfile), wait = TRUE,
		intern = FALSE, ignore.stderr = TRUE), silent = TRUE)
	if (inherits(res, "try-error") || !length(res)) return(character(0))
	if (res > 0) return(character(0)) # User cancelled input
	res <- readLines(tfile)
	return(res)	
}

## Linux/Unix version
.unixDlgDir <- function (default = getwd(), title = "")
{
    ## zenity must be installed on this machine!
    if (Sys.which("zenity") == "") return(NULL)
    ## Avoid displaying warning message in case user clicks on Cancel
    owarn <- getOption("warn")
    on.exit(options(warn = owarn))
    options(warn = -1)
    ## Use zenity to display the directory selection
	## There is no message area here, but one can set the title
	if (title == "") {
		title <- "Choose a directory" # Default title
	} #else {
	#	## Determine if the title is multiline...
	#	if (regexpr("\n", title) > 0) {
	#		## Try to use a notification instead
	#		if (Sys.which("notify-send") != "") {
	#			system(paste("notify-send --category=\"R\"",
	#				" \"R message\" \"", title, "\"", sep = ""), wait = FALSE)
	#			title <- "Choose folder"			
	#		} # Else the wole title cannot be displayed!!
	#	}
	#}
    msg <- paste("zenity --file-selection --title=\"", title,
	"\" --directory --filename=\"", default, "\"", sep = "")
    res <- system(msg, intern = TRUE)
    return(res)	
}
