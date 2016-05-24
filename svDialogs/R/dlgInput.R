## Define the S3 method
dlgInput <- function (message = "Enter a value",
default = "", ..., gui = .GUI)
{
	if (!gui$startUI("dlgInput", call = match.call(), default = default,
		msg = "Displaying a modal input dialog box",
		msg.no.ask = "A modal input dialog box was by-passed"))
		return(invisible(gui))
	
	## Check and rework main arguments and place them in gui$args
	if (!length(message)) message <- "Enter a value"
	message <- paste(as.character(message), collapse = "\n")
    if (is.null(default)) default <- "" else default <- as.character(default)[1]
	gui$setUI(args = list(message = message, default = default))
	
	## ... and dispatch to the method
	UseMethod("dlgInput", gui)
}

## Used to break the chain of NextMethod(), searching for a usable method
## in the current context
dlgInput.gui <- function (message = "Enter a value",
default = "", ..., gui = .GUI)
{
	msg <- paste("No workable method available to display an input dialog box using:",
		paste(guiWidgets(gui), collapse = ", "))
	gui$setUI(status = "error", msg = msg, widgets = "none")
	stop(msg)
}

## The pure textual version used a fallback in case no GUI could be used
dlgInput.textCLI <- function (message = "Enter a value",
default = "", ..., gui = .GUI)
{
	gui$setUI(widgets = "textCLI")
	## Ask for the input with readline()
	res <- readline(paste(gui$args$message,
		" [", gui$args$default, "]: ", sep = ""))
	if (res == "") res <- gui$args$default
	gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}

## The native version of the input box
dlgInput.nativeGUI <- function (message = "Enter a value",
default = "", ..., gui = .GUI)
{
	gui$setUI(widgets = "nativeGUI")
	## A simple text input box using native window
	## Return either a string, or character(0) if 'Cancel' clicked
	res <- switch(Sys.info()["sysname"],
		Windows = .winDlgInput(gui$args$message, gui$args$default),
		Darwin = .macDlgInput(gui$args$message, gui$args$default),
		.unixDlgInput(gui$args$message, gui$args$default)
	)
	
	## Do we need to further dispatch?
	if (is.null(res)) NextMethod("dlgInput", gui) else {
		gui$setUI(res = res, status = NULL)
		return(invisible(gui))
	}
}

## Windows version
.winDlgInput <- function (message, default)
{
	res <- winDialogString(message = message, default = default)
    if (is.null(res)) res <- character(0)
	return(res)
}

## Mac OS X version
.macDlgInput <- function (message, default)
{
    ## Display a modal message with native Mac dialog box
	#if (.Platform$GUI == "AQUA") app <- "(name of application \"R\")" else
	## This works from Mac OS X 10.5 Leopard:
	if (.Platform$GUI == "AQUA") app <- "(name of application id \"Rgui\")" else
		if (.isJGR()) app <- "\"JGR\"" else app <- "\"Terminal\""
	## Avoid displaying warning message when the user clicks on 'Cancel'
	owarn <- getOption("warn")
	on.exit(options(warn = owarn))
	options(warn = -1)
	cmd <- paste("-e 'tell application ", app,
		" to display dialog \"", message, "\" default answer \"", default,
		"\" with title \"Question\" buttons {\"Cancel\",\"OK\"} cancel button 1",
		" default button 2'", sep = "")
	## FI cannot use system(intern = TRUE) with this in R.app/R64.app
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
	res <- sub("^text returned:", "", res)
	res <- sub(", button returned:.*$", "", res)
	res <- paste(res, collapse = " ")
	return(res)
}

## Linux/Unix version
.unixDlgInput <- function (message, default) 
{
    ## zenity must be installed on this machine!
    if (Sys.which("zenity") == "") return(NULL)
    ## Avoid displaying warning message in case user clicks on Cancel
    owarn <- getOption("warn")
    on.exit(options(warn = owarn))
    options(warn = -1)
    ## Use zenity to display the prompt box
    msg <- paste("zenity --entry --title=\"Question\" --text=\"", message,
	"\" --entry-text=\"", default, "\"", sep = "")
    res <- system(msg, intern = TRUE)
    return(res)
}
