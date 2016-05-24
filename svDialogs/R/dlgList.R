## Define the S3 method
dlgList <- function (choices, preselect = NULL, multiple = FALSE, title = NULL,
..., gui = .GUI)
{
	## Check arguments
	choices <- as.character(choices)
	if (!length(choices)) return(character(0)) # Nothing to select	
	preselect <- as.character(preselect)
	preselect <- preselect[preselect %in% choices]
	if (!length(preselect)) preselect <- choices[1] # Select first item by default
	
	## Start a GUI action... or by-pass it
	if (!gui$startUI("dlgList", call = match.call(), default = preselect,
		msg = "Displaying a modal list dialog box",
		msg.no.ask = "A modal list dialog box was by-passed"))
		return(invisible(gui))
	
	## Further argument checking
	multiple <- isTRUE(as.logical(multiple))
	if (!length(title)) title <- NULL else title <- as.character(title)[1]
	gui$setUI(args = list(choices = choices, preselect = preselect,
		multiple = multiple, title = title))
	
	## ... and dispatch to the method
	UseMethod("dlgList", gui)
}

## Used to break the chain of NextMethod(), searching for a usable method
## in the current context
dlgList.gui <- function (choices, preselect = NULL, multiple = FALSE,
title = NULL, ..., gui = .GUI)
{
	msg <- paste("No workable method available to display a list dialog box using:",
		paste(guiWidgets(gui), collapse = ", "))
	gui$setUI(status = "error", msg = msg, widgets = "none")
	stop(msg)
}

## The pure textual version used a fallback in case no GUI could be used
dlgList.textCLI <- function (choices, preselect = NULL, multiple = FALSE,
title = NULL, ..., gui = .GUI)
{
	gui$setUI(widgets = "textCLI")
	## Ask a selection in a textual menu
	choices <- gui$args$choices
	multiple <- gui$args$multiple
	res <- select.list(choices = choices, preselect = gui$args$preselect,
		multiple = multiple, title = gui$args$title, graphics = FALSE)
    ## When multiple is FALSE and user cancelled, returns "" instead of
	## character(0) => change this for consistency
	if (!multiple && res == "" && !"" %in% choices)
		res <- character(0)
	gui$setUI(res = res, status = NULL)
	return(invisible(gui))
}

## The native version of the list box
dlgList.nativeGUI <- function (choices, preselect = NULL, multiple = FALSE,
title = NULL, ..., gui = .GUI)
{
	gui$setUI(widgets = "nativeGUI")
	## This is a simple 'select in the list' dialog box
    ## It follows the syntax of the select.list() function
	res <- switch(Sys.info()["sysname"],
		Windows = .winDlgList(gui$args$choices, gui$args$preselect,
			gui$args$multiple, gui$args$title),
		Darwin = .macDlgList(gui$args$choices, gui$args$preselect,
			gui$args$multiple, gui$args$title),
		.unixDlgList(gui$args$choices, gui$args$preselect,
			gui$args$multiple, gui$args$title)
	)
	
	## Do we need to further dispatch?
	if (is.null(res)) NextMethod("dlgList", gui) else {
		gui$setUI(res = res, status = NULL)
		return(invisible(gui))
	}
}

## Windows version
## select.list() does not have exactly the same behaviour, or native look&feel
## on all R programs => redefine it here to make it more homogeneous
.winDlgList <- function (choices, preselect = NULL, multiple = FALSE,
title = NULL)
{
	## To get the same behaviour as under Mac and Linux, that is, to return
	## character(0) when the user clicks on 'Cancel'
	res <- select.list(choices = choices, preselect = preselect,
		multiple = multiple, title = title, graphics = TRUE)
	if (length(res) == 1 && res == "") return(character(0)) else return(res)
}

## Mac OS X version
.macDlgList <- function (choices, preselect = NULL, multiple = FALSE,
title = NULL)
{
	if (.Platform$GUI == "AQUA") { # Use the R(64).app list box
		## Same as select.list(), but force graphics to TRUE
		## and always do preselection (first element) to match behaviour
		## under Windows
		if (is.null(preselect) || !any(preselect %in% choices))
			preselect <- choices[1]
		return(select.list(choices = choices, preselect = preselect,
			multiple = multiple, title = title, graphics = TRUE))
	} else if (.isJGR()) { # This seems to be JGR
		app <- "JGR"
	} else app <- "Terminal" # Probably run from terminal
	## Use osascript to display the list box
	## Make sure to keep only first preselection if !multiple
	if (!multiple) preselect <- preselect[1]
	## Format preselect into a single character string
	sel <- paste('"', preselect, '  "', sep = "", collapse = ",")
	## Format choices in a single string
	items <- paste('"', choices, '  "', sep = "", collapse = ",")
	## Default title
	if (is.null(title)) if (multiple) title <- "Select one or more" else
		title <- "Select one"
	## Avoid displaying warning message when the user clicks on 'Cancel'
	owarn <- getOption("warn")
	on.exit(options(warn = owarn))
	options(warn = -1)
	cmd <- paste("-e 'tell application \"", app, "\" to choose from list {",
		items, "} with title \"Make your selection\" with prompt \"", title,
		"\" multiple selections allowed ", multiple, " default items {",
		sel, "}'", sep = "")
	#res <- system2("osascript", cmd, stdout = TRUE, stderr = TRUE, wait = TRUE)
	res <- system(paste("osascript", cmd), intern = TRUE, wait = TRUE)
	if (res == "false") {
		return(character(0))
	} else {
		res <- unlist(strsplit(sub("  $", "", res), "  , ", fixed = TRUE))
		return(res)
	}
}

## Linux/Unix version
.unixDlgList <- function (choices, preselect = NULL, multiple = FALSE,
title = NULL)
{
    ## We don't use the ugly (on Linux) Tk version tk_select.list()
    ## In zenity, the normal list mode do not allow for preselections
    ## => switch to --checklist (multiple) or --radiolist (single) in this case
    ## Make sure that we have at least one preselection
    if (multiple) {
		kind <- "--checklist --column=\"Pick\" --column=\"Item\" --multiple"
	} else {
		kind <- "--radiolist --column=\"Pick\" --column=\"Item\""
	}
	## Only one item is preselected if multiple is FALSE (keep first one)
    if (!multiple) preselect <- preselect[1]
	## Create a string with TRUE/FALSE item alternated
    sel <- choices %in% preselect
    items <- paste(sel, ' "', choices, '"', sep = "", collapse = " ")
    ## Default title
    if (is.null(title)) if (multiple) title <- "Select one or more" else
	title <- "Select one"
    ## Avoid warning when user clicks on 'Cancel'
    owarn <- getOption("warn")
    on.exit(options(warn = owarn))
    options(warn = -1)
    ## Construct the command to send to zenity
    cmd <- paste("zenity --list --text=\"", title, "\" ", kind,
	" --hide-header --title=\"Make your choice\" --separator=\"@@@\" --height=",
	80 + 25 * length(choices), " ", items, sep = "")
    res <- system(cmd, intern = TRUE)
    res <- unlist(strsplit(res, "@@@", fixed = TRUE))
    if (is.null(res)) return(character(0)) else return(res)
}
